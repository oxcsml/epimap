// Constants
const TOPOJSON_PATH = "uk_lad_boundaries.json";
const RT_PATH = "Rt.csv";
const SITE_DATA_PATH = "site_data.csv";
const CASE_PROJECTION_PATH = "Cproj.csv";
const NHS_SCOTLAND_MAP = "nhs_scotland_health_boards.csv";
const ENGLAND_META_AREA_MAP = "england_meta_areas.csv";

// Set up dimensions for map
var map_svg = d3.select("#map"),
    width = +map_svg.attr("width"),
    height = +map_svg.attr("height");

var g = map_svg.append("g");
var barHeight = 20;
var barWidth = Math.floor(height / 3);
var margin = ({top: 20, right: 40, bottom: 30, left: 40})

// Set up dimensions for chart
var chart_margin = {top: 30, right: 30, bottom: 30, left: 30};
var chart_svg = d3.select("#chart"),
    chart_width = +chart_svg.attr("width") - chart_margin.left - chart_margin.right,
    chart_height = +chart_svg.attr("height") - chart_margin.top - chart_margin.bottom;
var chart_g = chart_svg.append("g")
    .attr("transform", "translate(" + chart_margin.left + "," + chart_margin.top + ")");

var actualChartLine = chart_svg.append("path")
    .attr("class", "actual-cases-line")

var smoothedChartLine = chart_svg.append("path")
    .attr("class", "smoothed-cases-line")

var projectedChartLine = chart_svg.append("path")
    .attr("class", "projected-cases-median-line");

var projectedArea = chart_svg.append("path")
  .attr("class", "projected-cases-area");

var casesLast7Info = d3.select("#cases-last7-info");
var casesLast7PerInfo = d3.select("#cases-last7-per-info");
var casesTotalInfo = d3.select("#cases-total-info");
var rtInfo = d3.select("#rt-info");
var caseProjInfo = d3.select("#case-proj-info");

// Add the X Axis
var chart_x_axis = chart_svg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chart_height})`);

// Add the Y Axis
var chart_y_axis = chart_svg.append("g")
    .attr("class", "chart-y-axis")
    .attr("transform", `translate(${chart_margin.left},0)`);

// Add a title
var chart_title = chart_svg.append("text")
    .attr("x", (chart_width / 2))             
    .attr("y", (chart_margin.top / 2))
    .attr("text-anchor", "middle")  
    .style("font-size", "16px");  

// Map and projection
var projection = d3.geoMercator()
        .center([-3.5, 54])
        .scale(3250)
        .translate([width/2, height/2]);
var path = d3.geoPath().projection(projection);

// Zooming
var zoomIn = map_svg.append("g").append("text")
    .attr("x", width-barHeight/2)
    .attr("y", Math.floor(height / 3) + margin.top + 70)
    .attr("width", 20)
    .attr("height", 20)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_in")
    .style("cursor", "pointer")
    .text("+");

var zoomIn = map_svg.append("g").append("text")
    .attr("x", width-barHeight/2)
    .attr("y", Math.floor(height / 3) + margin.top + 90)
    .attr("width", 20)
    .attr("height", 20)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_out")
    .style("cursor", "pointer")
    .text("-");

let zoom = d3.zoom()
    .on("zoom", ()=>g.selectAll("path").attr("transform", d3.event.transform));

map_svg.call(zoom)
    .on("dblclick.zoom", null);

d3.select("#zoom_in").on("click", function() {
    map_svg.transition().duration(500).call(zoom.scaleBy, 2);
});
d3.select("#zoom_out").on("click", function() {
    map_svg.transition().duration(500).call(zoom.scaleBy, 0.5);
});

// Data containers
var rtData = d3.map();
var caseTimeseries = d3.map();
var caseProjTimeseries = d3.map();
var nextWeekCaseProj = d3.map();
var caseHistory = d3.map();
var groupedAreaMap = d3.map();
var groupedAreaConstituents = d3.map();

// Tooltip container
var tooltip_div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);
var tooltip_header = tooltip_div.append("span")
    .attr("class", "header");
tooltip_div.append("br");
var tooltip_info1 = tooltip_div.append("span")
    .attr("class", "info-row1");
tooltip_div.append("br");
var tooltip_info2 = tooltip_div.append("span")
    .attr("class", "info-row2");
tooltip_div.append("br");
var tooltip_info3 = tooltip_div.append("span")
    .attr("class", "info-row3");

// Load external data
const loadCases = d3.csv(SITE_DATA_PATH).then(data=>{
    data.forEach(d=>{
        if (!caseTimeseries.has(d.area)) {
            caseTimeseries.set(d.area, []);
        }
        d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
        d.cases_new = +d.cases_new
        d.cases_new_smoothed = +d.cases_new_smoothed
        caseTimeseries.get(d.area).push(d);
    });
}).then(() => {
    caseTimeseries.each((cases, area) => {
        var casesLast7Day = cases.slice(1).slice(-7).map(c=>c.cases_new).reduce((a,b)=>a+b);
        var casesTotal = cases.map(c=>c.cases_new).reduce((a,b)=>a+b);
        caseHistory.set(area, {
            casesLast7Day: casesLast7Day,
            casesTotal: casesTotal
        });
    });
});

const urlParams = new URLSearchParams(window.location.search);
const rt_path = urlParams.get('rt') || RT_PATH;
const case_projection_path = urlParams.get('cproj') || CASE_PROJECTION_PATH;

const loadRt = d3.csv(rt_path).then(data => data.forEach(d => rtData.set(d.area, 
    {
        Rtlower: +d.Rt_lower,
        Rtmedian: +d.Rt_median,
        Rtupper: +d.Rt_upper,
    })));

const loadCaseProjections = d3.csv(case_projection_path).then(data => data.forEach(d => {
    if (!caseProjTimeseries.has(d.area)) {
        caseProjTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    d.C_lower = +d.C_lower;
    d.C_median = +d.C_median;
    d.C_upper = +d.C_upper;

    caseProjTimeseries.get(d.area).push(d);

})).then(() => {
    caseProjTimeseries.each((projections, area) => {
        var caseProjLower = 0, caseProjMedian = 0, caseProjUpper = 0;
        for (var i = 0; i < 7; i++) {
            caseProjLower += projections[i].C_lower;
            caseProjMedian += projections[i].C_median;
            caseProjUpper += projections[i].C_upper;
        }
        nextWeekCaseProj.set(area, {
            caseProjLower: Math.round(caseProjLower),
            caseProjMedian: Math.round(caseProjMedian),
            caseProjUpper: Math.round(caseProjUpper)
        });
    });
});

const loadNHSScotland = d3.csv(NHS_SCOTLAND_MAP).then(data => data.forEach(d => {
    const groupedArea = d["NHS Scotland Health Board"];
    groupedAreaMap.set(d.area, groupedArea);
    if (!groupedAreaConstituents.has(groupedArea)) {
        groupedAreaConstituents.set(groupedArea, []);
    }
    groupedAreaConstituents.get(groupedArea).push(d.area);
}));

const loadEnglandMetaAreas = d3.csv(ENGLAND_META_AREA_MAP).then(data => data.forEach(d => {
    const groupedArea = d["Meta area"];
    groupedAreaMap.set(d.area, groupedArea);
    if (!groupedAreaConstituents.has(groupedArea)) {
        groupedAreaConstituents.set(groupedArea, []);
    }
    groupedAreaConstituents.get(groupedArea).push(d.area);
}));

Promise.all([
    d3.json(TOPOJSON_PATH),
    loadRt,
    loadCaseProjections,
    loadCases,
    loadNHSScotland,
    loadEnglandMetaAreas
]).then(ready).catch(e=>{console.log("ERROR", e); throw e;});

var colorDomain = [0.5, 1.0, 2.0];

function getRtForArea(area) {
    var rt = rtData.get(area);
    var rtMedian = rt ? +rt.Rtmedian.toFixed(2) : "?";
    var rtUpper = rt ? +rt.Rtupper.toFixed(2) : "?";
    var rtLower = rt ? +rt.Rtlower.toFixed(2) : "?";
    return `${rtMedian} [${rtLower} - ${rtUpper}]`;
}

function getCaseProjForArea(area) {
    if (!nextWeekCaseProj.has(area)) {
        return "Unknown";
    }

    var projection = nextWeekCaseProj.get(area);

    var cprojmedian = projection.caseProjMedian;
    var cprojlower = projection.caseProjLower;
    var cprojupper = projection.caseProjUpper;

    return `${cprojmedian} [${cprojlower} - ${cprojupper}]`;
}

function getCaseHistoryForArea(area) {
    if (!caseHistory.has(area)) {
        return {
            casesLast7Day: "Unknown",
            casesTotal: "Unknown"
        };
    }
    return caseHistory.get(area);
}

// Handle data loaded
function ready(data) {
    var topo = data[0];

    console.log("Drawing map");

    // 1 is yellow, below 1 is green, above is red
    var rtColorScale = d3.scaleDiverging(t => d3.interpolateRdYlGn(1-t) )
        .domain(colorDomain);

    minCases = 1;
    maxCases = 50; // d3.max(nextWeekCaseProj.values().map(r=>r.caseProjMedian));
    console.log("minCases:", minCases, "maxCases:", maxCases);
    const logScale = d3.scaleLog().domain([minCases, maxCases]);
    const caseColorScale = d3.scaleSequential(v => d3.interpolateOrRd(logScale(v)));

    console.log("logScale.ticks:", logScale.ticks());

    var rtAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([colorDomain[0], colorDomain[2]]);

    var caseAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([minCases, maxCases]);

    var rtAxisFn = () => d3.axisBottom(rtAxisScale)
        .tickValues(rtColorScale.ticks())
        .tickFormat(rtColorScale.tickFormat())
        .tickSize(-barHeight);

    var caseAxisFn = () => d3.axisBottom(caseAxisScale)
        .tickValues(logScale.ticks(2))
        .tickFormat(d=>d)
        .tickSize(-barHeight);

    var axisBottom = map_svg.append("g")
        .attr("class", `x-axis`)
        .attr("transform", `translate(${width-barHeight},${margin.top}) rotate(90)`)
        .call(rtAxisFn)
        .selectAll("text")
        .attr("transform", "translate(-5, 15) rotate(-90)");

    rtFillFn = d => {  // Fill based on value of Rt
        var rt = rtData.get(d.properties.lad20nm);
        if (!rt) {
            return "#ccc";
        }
        return rtColorScale(rt.Rtmedian);
    }

    caseFillFn = d => { // Fill based on value of case projection
        var caseProj = nextWeekCaseProj.get(d.properties.lad20nm);
        if (!caseProj) {
            return "#ccc";
        }
        return caseColorScale(caseProj.caseProjMedian);
    }

    // Draw the map
    var map = g.selectAll("path")
        .data(topojson.feature(topo, topo.objects.Local_Authority_Districts__May_2020__Boundaries_UK_BFC).features)
        .enter().append("path")
        .attr("fill", rtFillFn)
        .style("fill-opacity", 1)
        .on("mouseover", function(d) {  // Add Tooltip on hover
            tooltip_div.transition()
              .duration(200)
              .style("opacity", .9);
            
            console.log(d.properties);

            tooltip_header.text(d.properties.lad20nm);
            tooltip_info1.text(`Last 7 days cases: ${getCaseHistoryForArea(d.properties.lad20nm).casesLast7Day}`);
            tooltip_info2.text(`Rt: ${getRtForArea(d.properties.lad20nm)}`);
            tooltip_info3.text(`Projected Cases: ${getCaseProjForArea(d.properties.lad20nm)}`);

            tooltip_div
              .style("left", (d3.event.pageX + 20) + "px")             
              .style("top", (d3.event.pageY - 28) + "px");
            d3.select(this).style("fill-opacity", 0.5);
        })
        .on("mouseout", function(d) {
            tooltip_div.style("opacity", 0);
            d3.select(this).style("fill-opacity", 1)
        })
        .on("click", d => selectArea(d.properties.lad20nm))
        .attr("d", path)
        .attr("class", "feature")

    // g.append("path")
    //     .datum(topojson.mesh(topo, topo.objects.Local_Authority_Districts__May_2020__Boundaries_UK_BFC, (a, b) => a !== b ))
    //     .attr("class", "mesh")
    //     .attr("d", path);

    // Draw the color scale
    const defs = map_svg.append("defs");
  
    const rtGradient = defs.append("linearGradient")
        .attr("id", "rt-gradient");

    const caseGradient = defs.append("linearGradient")
        .attr("id", "case-gradient");

    rtGradient.selectAll("stop")
        .data(rtColorScale.ticks().map((t, i, n) => ({ offset: `${100*i/n.length}%`, color: rtColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    caseGradient.selectAll("stop")
        .data(logScale.ticks().map((t, i, n) => ({ offset: `${100*i/n.length}%`, color: caseColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    var legend = map_svg.append('g')
        .attr("transform", `translate(${width},${barHeight}) rotate(90)`)
        .append("rect")
        .attr('transform', `translate(${margin.left}, 0)`)
        .attr("width", barWidth)        
        .attr("height", barHeight)
        .style("fill", "url(#rt-gradient)");

    legendText = map_svg.append("text")
        .attr("x", width-barHeight/2)       
        .attr("y", margin.top + 30)
        .attr("text-anchor", "middle")
        .style("font-size", "12px")  
        .text("Rt");

    // Add Rt vs case projection selection
    var showRt = map_svg.append("text")
        .attr("x", margin.left)             
        .attr("y", margin.top + 10)
        .style("font-size", "16px")  
        .style("cursor", "pointer")
        .attr("class", "active")
        .text("Rt");

    var showCases = map_svg.append("text")
        .attr("x", margin.left)             
        .attr("y", margin.top + 30)
        .style("font-size", "16px")  
        .style("cursor", "pointer")
        .text("Case Projections");

    showRt.on("click", () => {
        if (showCases.classed("active")) {
            map.attr("fill", rtFillFn);
            legend.style("fill", "url(#rt-gradient)");
            axisBottom.call(rtAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("Rt");
        }
        showRt.attr("class", "active");
        showCases.attr("class", "");
    });

    showCases.on("click", () => {
        if (showRt.classed("active")) {
            map.attr("fill", caseFillFn);
            legend.style("fill", "url(#case-gradient)");
            axisBottom.call(caseAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("Cases");
        }
        showCases.attr("class", "active");
        showRt.attr("class", "");
    });
    
}

function selectArea(selectedArea) {
    let area = selectedArea;
    if (groupedAreaMap.has(selectedArea)) {
        area = groupedAreaMap.get(selectedArea);
        const otherAreas = groupedAreaConstituents.get(area).join(", ");
        d3.select("#sub-heading").text(`Data shown is for the larger reporting area, ${area}, which contains ${otherAreas}`);
        d3.select("#cases-title").text(`Cases for ${area} (including ${selectedArea})`);
        d3.select("#estimates-title").text(`Estimates for ${area} (including ${selectedArea})`);
    }
    else {
        d3.select("#sub-heading").text("");
        d3.select("#cases-title").text(`Cases`);
        d3.select("#estimates-title").text(`Estimates`);
    }

    d3.select("#data-heading").text(selectedArea);
    

    var chartData = caseTimeseries.get(area);
    if (!chartData) {
        console.log("ERROR: No chart data found for area ", area);
        return;
    }
    var projectionData = caseProjTimeseries.get(area);
    if (!projectionData) {
        console.log("ERROR: No projection data found for area ", area);
        return;
    }
    console.log("Case projection data for " + area, projectionData);

    var xDomain = d3.extent([...chartData.map(c=>c.Date), ...projectionData.map(p=>p.Date)]);
    var yDomain = [0, 52]; //d3.max([...chartData.map(c=>c.cases_new), ...projectionData.map(p=>p.C_median)])];

    console.log("X Domain:", xDomain);
    console.log("Y Domain:", yDomain);

    var x = d3.scaleTime()
        .domain(xDomain)
        .range([chart_margin.left, chart_width-chart_margin.right]);
    var y = d3.scaleLinear()
        .domain(yDomain)
        .range([chart_height, 0]);

    // Define the lines
    var actualCasesLine = d3.line()
        .x(function(d) { return x(d.Date); })
        .y(function(d) { return y(d.cases_new); });

    var smoothedCasesLine = d3.line()
        .x(function(d) { return x(d.Date); })
        .y(function(d) { return y(d.cases_new_smoothed); });

    var projectedCasesLine = d3.line()
        .x(function(d) { return x(d.Date); })
        .y(function(d) { return y(d.C_median); });

    var projectedCasesArea = d3.area()
        .x(function(d) { return x(d.Date); })
        .y0(function(d) { return y(d.C_lower); })
        .y1(function(d) { return y(d.C_upper); });

    actualChartLine
        .datum(chartData)
        .transition()
        .duration(500)
        .attr("d", actualCasesLine);

    smoothedChartLine
        .datum(chartData)
        .transition()
        .duration(500)
        .attr("d", smoothedCasesLine);

    projectedArea
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedCasesArea);

    projectedChartLine
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedCasesLine);

    var caseHistory = getCaseHistoryForArea(area);
    casesLast7Info.text(caseHistory.casesLast7Day);
    casesTotalInfo.text(caseHistory.casesTotal);
    rtInfo.text(getRtForArea(area));
    caseProjInfo.text(getCaseProjForArea(area));

    chart_x_axis.call(d3.axisBottom(x));
    chart_y_axis.call(d3.axisLeft(y).ticks(5));
    chart_title.text(`COVID-19 Cases for ${area}`);
    
}

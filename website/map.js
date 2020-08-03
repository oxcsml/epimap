// Constants
const TOPOJSON_PATH = "uk_boundaries.json";
const RT_PATH = "Rt.csv";
const SITE_DATA_PATH = "site_data.csv";
const CASE_PROJECTION_PATH = "Cproj.csv";

// Set up dimensions for map
var map_svg = d3.select("#map"),
    width = +map_svg.attr("width"),
    height = +map_svg.attr("height");

var g = map_svg.append("g");
var barHeight = 20;
var barWidth = Math.floor(height / 3);
var margin = ({top: 20, right: 40, bottom: 30, left: 40})

// Set up dimensions for chart
var chart_margin = {top: 30, right: 30, bottom: 30, left: 60};
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
    .attr("transform", "translate(0," + chart_height + ")");

// Add the Y Axis
var chart_y_axis = chart_svg.append("g")
    .attr("class", "chart-y-axis");

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

// Data containers
var rtData = d3.map();
var caseTimeseries = d3.map();
var caseProjTimeseries = d3.map();
var nextWeekCaseProj = d3.map();

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
var loadCases = d3.csv(SITE_DATA_PATH).then(data=>{
    data.forEach(d=>{
        if (!caseTimeseries.has(d.area)) {
            caseTimeseries.set(d.area, []);
        }
        d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
        d.cases_new = +d.cases_new
        d.cases_new_smoothed = +d.cases_new_smoothed
        caseTimeseries.get(d.area).push(d);
    });
})

var loadRt = d3.csv(RT_PATH).then(data => data.forEach(d => rtData.set(d.area, 
    {
        Rtlower: +d.Rt_lower,
        Rtmedian: +d.Rt_median,
        Rtupper: +d.Rt_upper,
    })));

var loadCaseProjections = d3.csv(CASE_PROJECTION_PATH).then(data => data.forEach(d => {
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
            caseProjLower: caseProjLower,
            caseProjMedian: caseProjMedian,
            caseProjUpper: caseProjUpper
        });
    });
});

Promise.all([
    d3.json(TOPOJSON_PATH),
    loadRt,
    loadCaseProjections,
    loadCases
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

    var cprojmedian = projection.caseProjLower.toFixed(2);
    var cprojlower = projection.caseProjMedian.toFixed(2);
    var cprojupper = projection.caseProjUpper.toFixed(2);

    return `${cprojmedian} [${cprojlower} - ${cprojupper}]`;
}

// Handle data loaded
function ready(data) {
    var topo = data[0];

    console.log("Drawing map");

    // 1 is yellow, below 1 is green, above is red
    var rtColorScale = d3.scaleDiverging(t => d3.interpolateRdYlGn(1-t) )
        .domain(colorDomain);

    minCases = 1;
    maxCases = 100; // d3.max(nextWeekCaseProj.values().map(r=>r.caseProjMedian));
    console.log("minCases:", minCases, "maxCases:", maxCases);
    const logScale = d3.scaleLog().domain([minCases, maxCases]);
    const caseColorScale = d3.scaleSequential(v => d3.interpolateOrRd(logScale(v)));

    console.log("Rt color scale domain:", rtColorScale.domain())

    // TODO: Replace legend color scale with cases scale when selected.
    var axisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([colorDomain[0], colorDomain[2]]);

    var axisBottom = g => g
        .attr("class", `x-axis`)
        .attr("transform", `translate(${width-barHeight},${margin.top}) rotate(90)`)
        .call(
            d3.axisBottom(axisScale)
                .tickValues(rtColorScale.ticks())
                .tickFormat(rtColorScale.tickFormat())
                .tickSize(-barHeight)
        )
        .selectAll("text")
        .attr("transform", "translate(-5, 15) rotate(-90)");

    rtFillFn = d => {  // Fill based on value of Rt
        var rt = rtData.get(d.properties.ctyua17nm);
        if (!rt) {
            return "#ccc";
        }
        return rtColorScale(rt.Rtmedian);
    }

    caseFillFn = d => { // Fill based on value of case projection
        var caseProj = nextWeekCaseProj.get(d.properties.ctyua17nm);
        if (!caseProj) {
            return "#ccc";
        }
        return caseColorScale(caseProj.caseProjMedian);
    }

    // Draw the map
    var map = g.selectAll("path")
        //.data(topojson.feature(topo, topo.objects.Counties_and_Unitary_Authorities__December_2016__Boundaries).features)
        .data(topojson.feature(topo, topo.objects.Counties_and_Unitary_Authorities__December_2017__Boundaries_UK).features)
        .enter().append("path")
        .attr("fill", rtFillFn)
        .style("fill-opacity", 1)
        .on("mouseover", function(d) {  // Add Tooltip on hover
            tooltip_div.transition()
              .duration(200)
              .style("opacity", .9);
            
            console.log(d.properties);

            tooltip_header.text(d.properties.ctyua17nm);
            tooltip_info1.text(`Last 7 days cases: TODO`);
            tooltip_info2.text(`Rt: ${getRtForArea(d.properties.ctyua17nm)}`);
            tooltip_info3.text(`Projected Cases: ${getCaseProjForArea(d.properties.ctyua17nm)}`);

            tooltip_div
              .style("left", (d3.event.pageX + 10) + "px")             
              .style("top", (d3.event.pageY - 28) + "px");
            d3.select(this).style("fill-opacity", 0.5);
        })
        .on("mouseout", function(d) {
            tooltip_div.style("opacity", 0);
            d3.select(this).style("fill-opacity", 1)
        })
        .on("click", d => selectArea(d.properties.ctyua17nm))
        .attr("d", path)
        .attr("class", "feature")

    g.append("path")
        .datum(topojson.mesh(topo, topo.objects.Counties_and_Unitary_Authorities__December_2017__Boundaries_UK, (a, b) => a !== b ))
        .attr("class", "mesh")
        .attr("d", path);

    // Draw the color scale
    const defs = map_svg.append("defs");
  
    const linearGradient = defs.append("linearGradient")
        .attr("id", "linear-gradient");

    linearGradient.selectAll("stop")
        .data(rtColorScale.ticks().map((t, i, n) => ({ offset: `${100*i/n.length}%`, color: rtColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    map_svg.append('g')
        .attr("transform", `translate(${width},${barHeight}) rotate(90)`)
        .append("rect")
        .attr('transform', `translate(${margin.left}, 0)`)
        .attr("width", barWidth)        
        .attr("height", barHeight)
        .style("fill", "url(#linear-gradient)");

    map_svg.append('g')
        .call(axisBottom);

    map_svg.append("text")
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
        .text("Rt")
        .on("click", function() {
            d3.select(this).attr("class", "active");
        });

    var showCases = map_svg.append("text")
        .attr("x", margin.left)             
        .attr("y", margin.top + 30)
        .style("font-size", "16px")  
        .style("cursor", "pointer")
        .text("Case Projections");

    showRt.on("click", () => {
        if (showCases.classed("active")) {
            map.attr("fill", rtFillFn);
            // TODO: Switch the color scale
        }
        showRt.attr("class", "active");
        showCases.attr("class", "");
    });

    showCases.on("click", () => {
        if (showRt.classed("active")) {
            map.attr("fill", caseFillFn);
        }
        showCases.attr("class", "active");
        showRt.attr("class", "");
    });
    
}

function selectArea(area) {
    d3.select("#data-heading").text(area);

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
    var yDomain = [0, 100]; //d3.max([...chartData.map(c=>c.cases_new), ...projectionData.map(p=>p.C_median)])];

    console.log("X Domain:", xDomain);
    console.log("Y Domain:", yDomain);

    var x = d3.scaleTime()
        .domain(xDomain)
        .range([0, chart_width]);
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

    rtInfo.text(getRtForArea(area));
    caseProjInfo.text(getCaseProjForArea(area));

    chart_x_axis.call(d3.axisBottom(x));
    chart_y_axis.call(d3.axisLeft(y).ticks(5));
    chart_title.text(`COVID-19 Cases for ${area}`);
    
}

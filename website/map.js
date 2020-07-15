// Constants
// const TOPOJSON_PATH = "https://raw.githubusercontent.com/rs-delve/Rmap/master/uk_boundaries.json?token=ABFV53BX3YW3FJDUK4BFBFS7C4RJU";
// const DATA_PATH = "https://raw.githubusercontent.com/rs-delve/Rmap/master/RtCproj.csv?token=ABFV53DIFS3BWCRSDDTZKFK7C4KQW"
// const SITE_DATA_PATH = 'https://raw.githubusercontent.com/rs-delve/Rmap/master/site_data.csv?token=ABFV53FFK2TGLFOINXD6CEK7DASKQ'

const TOPOJSON_PATH = "uk_boundaries.json";
const DATA_PATH = "RtCproj.csv";
const SITE_DATA_PATH = "site_data.csv";

// Set up dimensions for map
var map_svg = d3.select("#map"),
    width = +map_svg.attr("width"),
    height = +map_svg.attr("height");

var g = map_svg.append("g");
var barHeight = 30;
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

var casesLast7Info = d3.select("#cases-last7-info");
var casesLast7PerInfo = d3.select("#cases-last7-per-info");
var casesTotalInfo = d3.select("#cases-total-info");
var rtInfo = d3.select("#rt-info");

// Map and projection
var projection = d3.geoMercator()
        .center([-1.6,52.6])
        .scale(4000)
        .translate([width/2, height/2]);
var path = d3.geoPath().projection(projection);

// Data containers
var rtData = d3.map();
var caseTimeseries = d3.map();

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
var load_cases = d3.csv(SITE_DATA_PATH).then(data=>{
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

Promise.all([
    d3.json(TOPOJSON_PATH),
    d3.csv(DATA_PATH).then(data => data.forEach(d => rtData.set(d.area, 
        {
            Rtlower: +d.Rtlower,
            Rtmedian: +d.Rtmedian,
            Rtupper: +d.Rtupper,
            Cprojlower: +d.Cprojlower,
            Cprojmedian: +d.Cprojmedian,
            Cprojupper: +d.Cprojupper
        }))),
    load_cases
]).then(ready).catch(e=>{console.log("ERROR", e); throw e;});

var colorDomain = [0, 1, 4];

function getRtForArea(area) {
    var rt = rtData.get(area);
    var rtMedian = rt ? +rt.Rtmedian.toFixed(2) : "?";
    var rtUpper = rt ? +rt.Rtupper.toFixed(2) : "?";
    var rtLower = rt ? +rt.Rtlower.toFixed(2) : "?";
    return `${rtMedian} [${rtLower} - ${rtUpper}]`;
}

function getCaseProjForArea(area) {
    var rt = rtData.get(area);
    var Cprojmedian = rt ? +rt.Cprojmedian.toFixed(2) : "?";
    var Cprojlower = rt ? +rt.Cprojlower.toFixed(2) : "?";
    var Cprojupper = rt ? +rt.Cprojupper.toFixed(2) : "?";
    return `${Cprojmedian} [${Cprojlower} - ${Cprojupper}]`;
}

// Handle data loaded
function ready(data) {
    var topo = data[0];

    console.log("Drawing map");

    // 1 is yellow, below 1 is green, above is red
    var rtColorScale = d3.scaleDiverging(t => d3.interpolateRdYlGn(1-t) )
        .domain(colorDomain);

    minCases = 1;
    maxCases = d3.max(rtData.values().map(r=>r.Cprojmedian));
    console.log("minCases:", minCases, "maxCases:", maxCases);
    const logScale = d3.scaleLog().domain([minCases, maxCases]);
    const caseColorScale = d3.scaleSequential(v => d3.interpolateOrRd(logScale(v)));

    console.log("Rt color scale domain:", rtColorScale.domain())

    // TODO: Replace legend color scale with cases scale when selected.
    var axisScale = d3.scaleLinear()
        .range([margin.left, width - margin.right])
        .domain([colorDomain[0], colorDomain[2]]);

    var axisBottom = g => g
        .attr("class", `x-axis`)
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(
            d3.axisBottom(axisScale)
                .tickValues(rtColorScale.ticks())
                .tickFormat(rtColorScale.tickFormat())
                .tickSize(-barHeight)
        );

    rtFillFn = d => {  // Fill based on value of Rt
        var rt = rtData.get(d.properties.ctyua16nm);
        if (!rt) {
            return "#ccc";
        }
        return rtColorScale(rt.Rtmedian);
    }

    caseFillFn = d => { // Fill based on value of case projection
        var rt = rtData.get(d.properties.ctyua16nm);
        if (!rt) {
            return "#ccc";
        }
        return caseColorScale(rt.Cprojmedian);
    }

    // Draw the map
    var map = g.selectAll("path")
        .data(topojson.feature(topo, topo.objects.Counties_and_Unitary_Authorities__December_2016__Boundaries).features)
        .enter().append("path")
        .attr("fill", rtFillFn)
        .style("fill-opacity", 1)
        .on("mouseover", function(d) {  // Add Tooltip on hover
            tooltip_div.transition()
              .duration(200)
              .style("opacity", .9);
            
            tooltip_header.text(d.properties.ctyua16nm);
            tooltip_info1.text(`Last 7 days cases: TODO`);
            tooltip_info2.text(`Rt: ${getRtForArea(d.properties.ctyua16nm)}`);
            tooltip_info3.text(`Projected Cases: ${getCaseProjForArea(d.properties.ctyua16nm)}`);

            tooltip_div
              .style("left", (d3.event.pageX + 10) + "px")             
              .style("top", (d3.event.pageY - 28) + "px");
            d3.select(this).style("fill-opacity", 0.5);
        })
        .on("mouseout", function(d) {
            tooltip_div.style("opacity", 0);
            d3.select(this).style("fill-opacity", 1)
        })
        .on("click", d => selectArea(d.properties.ctyua16nm))
        .attr("d", path)
        .attr("class", "feature")

    g.append("path")
        .datum(topojson.mesh(topo, topo.objects.Counties_and_Unitary_Authorities__December_2016__Boundaries, (a, b) => a !== b ))
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
        .attr("transform", `translate(0,${height - margin.bottom - barHeight})`)
        .append("rect")
        .attr('transform', `translate(${margin.left}, 0)`)
        .attr("width", width - margin.right - margin.left)
        .attr("height", barHeight)
        .style("fill", "url(#linear-gradient)");
      
    map_svg.append('g')
        .call(axisBottom);

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

    chart_data = caseTimeseries.get(area);
    if (!chart_data) {
        console.log("ERROR: No chart data found for area ", area);
        return;
    }
    
    var x = d3.scaleTime()
        .domain(d3.extent(chart_data, function(d) { return d.Date; }))
        .range([0, chart_width]);
    var y = d3.scaleLinear()
        .domain([0, d3.max(chart_data, function(d) { return d.cases_new; })])
        .range([chart_height, 0]);

    // Define the lines
    var actualCasesLine = d3.line()
        .x(function(d) { return x(d.Date); })
        .y(function(d) { return y(d.cases_new); });

    var smoothedCasesLine = d3.line()
        .x(function(d) { return x(d.Date); })
        .y(function(d) { return y(d.cases_new_smoothed); });

    actualChartLine
        .datum(chart_data)
        .transition()
        .duration(500)
        .attr("d", actualCasesLine);

    smoothedChartLine
        .datum(chart_data)
        .transition()
        .duration(500)
        .attr("d", smoothedCasesLine);

    rtInfo.text(getRtForArea(area));

    // TODO: Move these out of the update function
    // Add the X Axis
    chart_svg.append("g")
        .attr("class", "chart-x-axis")
        .attr("transform", "translate(0," + chart_height + ")")
        .call(d3.axisBottom(x));

    // Add the Y Axis
    chart_svg.append("g")
        .attr("class", "chart-y-axis")
        .call(d3.axisLeft(y).ticks(5));

    // Add a title
    chart_svg.append("text")
        .attr("x", (chart_width / 2))             
        .attr("y", (chart_margin.top / 2))
        .attr("text-anchor", "middle")  
        .style("font-size", "16px")  
        .text(`COVID-19 Cases for ${area}`);
}
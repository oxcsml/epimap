// Constants
const TOPOJSON_PATH = "uk_lad_boundaries.json";
const MAP_PATH = 'default';
const RT_PATH = "Rt.csv";
const SITE_DATA_PATH = "site_data.csv";
const CASE_PROJECTION_PATH = "Cproj.csv";
const CASE_PREDICTION_PATH = "Cpred.csv";
const PEXCEED_PATH = "Pexceed.csv";
const NHS_SCOTLAND_MAP = "nhs_scotland_health_boards.csv";
const ENGLAND_META_AREA_MAP = "england_meta_areas.csv";
const METADATA_PATH = "metadata.csv";

// Set up dimensions for map
const MAX_CASES = 50;

const map_svg = d3.select("#map"),
    width = +map_svg.attr("width"),
    height = +map_svg.attr("height");

const barHeight = 20;
const sliderWidth = 300;
const barWidth = Math.floor(height / 3);
const margin = ({ top: 20, right: 40, bottom: 30, left: 40 });
const g = map_svg.append("g");

const sliderLeft = Math.round((width - margin.left - margin.right - sliderWidth) / 2);
const sliderSvg = d3.select("#slider-svg").style("display", "block");
const sliderG = sliderSvg.append("g")
    .attr('transform', `translate(${sliderLeft},10)`);
sliderSvg.append("text")
    .attr("x", sliderLeft-50)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .style("font-size", "13px")
    .text("Rt Date");
const sliderValueLabel = sliderSvg.append("text")
    .attr("x", sliderLeft +sliderWidth + 60)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .style("font-size", "13px");

// Set up dimensions for chart
const chartMargin = { top: 30, right: 30, bottom: 30, left: 30 };
const caseChartSvg = d3.select("#chart");
// Note: Assuming the two charts are the same size! 
const chartWidth = +caseChartSvg.attr("width") - chartMargin.left - chartMargin.right;
const chartHeight = +caseChartSvg.attr("height") - chartMargin.top - chartMargin.bottom;

const chartG = caseChartSvg.append("g")
    .attr("transform", "translate(" + chartMargin.left + "," + chartMargin.top + ")");

const actualChartLine = caseChartSvg.append("path")
    .attr("class", "actual-cases-line")

const smoothedChartLine = caseChartSvg.append("path")
    .attr("class", "smoothed-cases-line")

const predictedChartLine = caseChartSvg.append("path")
    .attr("class", "predicted-cases-median-line");

const predictedArea = caseChartSvg.append("path")
    .attr("class", "predicted-cases-area");

const projectedChartLine = caseChartSvg.append("path")
    .attr("class", "projected-cases-median-line");

const projectedArea = caseChartSvg.append("path")
    .attr("class", "projected-cases-area");

const rtChartSvg = d3.select("#rt-chart");

const rtChartMedianLine = rtChartSvg.append("path")
    .attr("class", "rt-line");

const rtChartInnerArea = rtChartSvg.append("path")
    .attr("class", "rt-inner-area");

const rtChartOuterArea = rtChartSvg.append("path")
    .attr("class", "rt-outer-area");

const rtHorizontalLine = rtChartSvg.append("g")
    .attr("class", "rt-horizontal-line");

// Add the X Axes
const caseChartXAxis = caseChartSvg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chartHeight})`);

const rtChartXAxis = rtChartSvg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chartHeight})`);

// Add the Y Axis
const caseChartYAxis = caseChartSvg.append("g")
    .attr("class", "chart-y-axis")
    .attr("transform", `translate(${chartMargin.left},0)`);

const rtChartYAxis = rtChartSvg.append("g")
    .attr("class", "chart-y-axis")
    .attr("transform", `translate(${chartMargin.left},0)`);

// Add a title
const caseChartTitle = caseChartSvg.append("text")
    .attr("x", (chartWidth / 2))
    .attr("y", (chartMargin.top / 2))
    .attr("text-anchor", "middle")
    .style("font-size", "16px");

const rtChartTitle = rtChartSvg.append("text")
    .attr("x", (chartWidth / 2))
    .attr("y", (chartMargin.top / 2))
    .attr("text-anchor", "middle")
    .style("font-size", "16px");

const casesLast7Info = d3.select("#cases-last7-info");
const casesLast7PerInfo = d3.select("#cases-last7-per-info");
const casesTotalInfo = d3.select("#cases-total-info");
const rtInfo = d3.select("#rt-info");
const caseProjInfo = d3.select("#case-proj-info");
const caseProjPer100kInfo = d3.select("#case-proj-per100k-info");

// Map and projection
const projection = d3.geoMercator()
    .center([-3.5, 54])
    .scale(3250)
    .translate([width / 2, height / 2]);
const path = d3.geoPath().projection(projection);

// Zooming
const zoomIn = map_svg.append("g").append("text")
    .attr("x", width - barHeight / 2)
    .attr("y", Math.floor(height / 3) + margin.top + 70)
    .attr("width", 20)
    .attr("height", 20)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_in")
    .style("cursor", "pointer")
    .text("+");

const zoomOut = map_svg.append("g").append("text")
    .attr("x", width - barHeight / 2)
    .attr("y", Math.floor(height / 3) + margin.top + 90)
    .attr("width", 20)
    .attr("height", 20)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_out")
    .style("cursor", "pointer")
    .text("-");

let zoom = d3.zoom()
    .on("zoom", () => g.selectAll("path").attr("transform", d3.event.transform));

map_svg.call(zoom)
    .on("dblclick.zoom", null);

d3.select("#zoom_in").on("click", function () {
    map_svg.transition().duration(500).call(zoom.scaleBy, 2);
});
d3.select("#zoom_out").on("click", function () {
    map_svg.transition().duration(500).call(zoom.scaleBy, 0.5);
});

const dateFormat = d3.timeFormat("%Y-%m-%d");

// Data containers
const rtData = d3.map();
const caseTimeseries = d3.map();
const caseProjTimeseries = d3.map();
const casePredTimeseries = d3.map();
const pexceedData = d3.map();
const nextWeekCaseProj = d3.map();
const nextWeekCaseProjPer100k = d3.map();
const caseHistory = d3.map();
const caseHistoryPer100k = d3.map();
const groupedAreaMap = d3.map();
const groupedAreaConstituents = d3.map();
const populations = d3.map();

// Tooltip container
const tooltip_div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);
const tooltip_header = tooltip_div.append("span")
    .attr("class", "header");
tooltip_div.append("br");
const tooltip_info1 = tooltip_div.append("span")
    .attr("class", "info-row1");
tooltip_div.append("br");
const tooltip_info2 = tooltip_div.append("span")
    .attr("class", "info-row2");
tooltip_div.append("br");
const tooltip_info3 = tooltip_div.append("span")
    .attr("class", "info-row3");
tooltip_div.append("br");
const tooltip_info4 = tooltip_div.append("span")
    .attr("class", "info-row4");

// Load external data
const loadCases = d3.csv(SITE_DATA_PATH).then(data => {
    data.forEach(d => {
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
        const casesLast7Day = cases.slice(1).slice(-7).map(c => c.cases_new).reduce((a, b) => a + b);
        const casesTotal = cases.map(c => c.cases_new).reduce((a, b) => a + b);
        caseHistory.set(area, {
            casesLast7Day: casesLast7Day,
            casesTotal: casesTotal
        });
    });
});

const urlParams = new URLSearchParams(window.location.search);
const map_path = urlParams.get("map") || MAP_PATH
const rt_path = map_path.concat("/", RT_PATH);
const case_projection_path = map_path.concat("/", CASE_PROJECTION_PATH);
const case_prediction_path = map_path.concat("/", CASE_PREDICTION_PATH);
const pexceed_path = map_path.concat("/", PEXCEED_PATH);

const loadRt = d3.csv(rt_path).then(data => {
    data.forEach(d => {
        if (!rtData.has(d.area)) {
            rtData.set(d.area, []);
        }
        const current = {
            Date: d3.timeParse("%Y-%m-%d")(d.Date),
            Rt025: +d.Rt_2_5,
            Rt10: +d.Rt_10,
            Rt20: +d.Rt_20,
            Rt25: +d.Rt_25,
            Rt30: +d.Rt_30,
            Rt40: +d.Rt_40,
            Rt50: +d.Rt_50,
            Rt60: +d.Rt_60,
            Rt70: +d.Rt_70,
            Rt75: +d.Rt_75,
            Rt80: +d.Rt_80,
            Rt90: +d.Rt_90,
            Rt975: +d.Rt_97_5
        };
        rtData.get(d.area).push(current);
    });
});

const loadPExceed = d3.csv(pexceed_path).then(data => {
    data.forEach(d => {
        if (!pexceedData.has(d.area)) {
            pexceedData.set(d.area, []);
        }

        const current = {
            Date: d3.timeParse("%Y-%m-%d")(d.Date),
            P10: +d.P_10
            // Ignoring other fields, P_08,P_09,P_11,P_12,P_15,P_20 for now
        };
        pexceedData.get(d.area).push(current);
    });
});

const loadCaseProjections = d3.csv(case_projection_path).then(data => data.forEach(d => {
    if (!caseProjTimeseries.has(d.area)) {
        caseProjTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.C_025) {
      d.C_lower = +d.C_025;
      d.C_median = +d.C_50;
      d.C_upper = +d.C_975;
    } else {
      d.C_lower = +d.C_lower;
      d.C_median = +d.C_median;
      d.C_upper = +d.C_upper;
    }  

    caseProjTimeseries.get(d.area).push(d);

})).then(() => {
    caseProjTimeseries.each((projections, area) => {
        let caseProjLower = 0, caseProjMedian = 0, caseProjUpper = 0;
        for (let i = 0; i < 7; i++) {
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

const loadCasePredictions = d3.csv(case_prediction_path).then(data => data.forEach(d => {
    if (!casePredTimeseries.has(d.area)) {
        casePredTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.C_025) {
      d.C_lower = +d.C_025;
      d.C_median = +d.C_50;
      d.C_upper = +d.C_975;
    } else {
      d.C_lower = +d.C_lower;
      d.C_median = +d.C_median;
      d.C_upper = +d.C_upper;
    }  

    casePredTimeseries.get(d.area).push(d);

}));

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

const loadMetadata = d3.csv(METADATA_PATH).then(data => data.forEach(d => {
    populations.set(d.AREA, +d.POPULATION);
}));

const getPopulation = (area) => {
    if (groupedAreaMap.has(area)) {
        const groupedArea = groupedAreaMap.get(area);
        const constituents = groupedAreaConstituents.get(groupedArea);

        return constituents
            .map(a => populations.get(a))
            .reduce((a, b) => a + b, 0);
    }
    return populations.get(area);
}

const casesAndMeta = Promise.all([
    loadCaseProjections,
    loadCasePredictions,
    loadCases,
    loadMetadata,
    loadNHSScotland,
    loadEnglandMetaAreas
]).then(() => {
    nextWeekCaseProj.each((caseProj, area) => {
        const pop = getPopulation(area) / 100000;
        nextWeekCaseProjPer100k.set(area, {
            caseProjLower: Math.round(caseProj.caseProjLower / pop),
            caseProjMedian: Math.round(caseProj.caseProjMedian / pop),
            caseProjUpper: Math.round(caseProj.caseProjUpper / pop)
        });
    });

    caseHistory.each((cases, area) => {
        const pop = getPopulation(area) / 100000;
        caseHistoryPer100k.set(area, {
            casesLast7Day: Math.round(cases.casesLast7Day / pop),
            casesTotal: Math.round(cases.casesTotal / pop)
        });
    });
});

Promise.all([
    d3.json(TOPOJSON_PATH),
    loadRt,
    loadPExceed,
    casesAndMeta
]).then(ready).catch(e => { console.log("ERROR", e); throw e; });

const colorDomain = [0.5, 1.0, 2.0];
const pExceedColorDomain = [0, 0.5, 1.0];

function getRtForArea(area) {
    const rtSeries = rtData.get(area);
    if (!rtSeries) {
        return "? [? - ?]";
    }
    const [lastRt] = rtSeries.slice(-1)
    const median = lastRt.Rt50.toFixed(2);
    const upper = lastRt.Rt90.toFixed(2);
    const lower = lastRt.Rt10.toFixed(2);
    return `${median} [${lower} - ${upper}]`;
}

function getCaseProjForArea(area) {
    if (!nextWeekCaseProj.has(area)) {
        return "Unknown";
    }

    const projection = nextWeekCaseProj.get(area);

    const cprojmedian = projection.caseProjMedian;
    const cprojlower = projection.caseProjLower;
    const cprojupper = projection.caseProjUpper;

    return `${cprojmedian} [${cprojlower} - ${cprojupper}]`;
}

function getPExceedForArea(area) {
    if (!pexceedData.has(area)) {
        return "Unknown";
    }
    const pExceedSeries = pexceedData.get(area);
    const [lastPExceed] = pExceedSeries.slice(-1)
    const pRtOver1 = lastPExceed.P10.toFixed(2);
    return `${pRtOver1}`;
}

function getCaseProjPer100kForArea(area) {
    if (!nextWeekCaseProj.has(area)) {
        return "Unknown";
    }

    const projection = nextWeekCaseProjPer100k.get(area);

    const cprojmedian = projection.caseProjMedian;
    const cprojlower = projection.caseProjLower;
    const cprojupper = projection.caseProjUpper;

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

function getCaseHistoryPer100kForArea(area) {
    if (!caseHistoryPer100k.has(area)) {
        return {
            casesLast7Day: "Unknown",
            casesTotal: "Unknown"
        };
    }
    return caseHistoryPer100k.get(area);
}

function areaNameToId(areaName) {
    return areaName.replace(" and ", "-").replace(" upon ", "-").replace(",", "").replace(" ", "-");
}

// Handle data loaded
function ready(data) {
    const topo = data[0];

    console.log("Drawing map");

    const rtColorScale = d3.scaleDiverging(t => d3.interpolateRdBu(1 - t))
        .domain(colorDomain);
    const pExceedColorScale = d3.scaleDiverging(t => d3.interpolateRdBu(1 - t))
        .domain(pExceedColorDomain);

    minCases = 1;
    maxCases = MAX_CASES; // d3.max(nextWeekCaseProjPer100k.values().map(r=>r.caseProjMedian));
    maxColorCases = d3.max(nextWeekCaseProjPer100k.values().map(r=>parseFloat(r.caseProjMedian)));
    const logScale = d3.scaleLog().domain([minCases, maxCases]);
    const caseLogScale = d3.scaleLog().domain([minCases, maxColorCases]).range([margin.left, margin.left + barWidth]);
    const caseColorScale = d3.scaleSequential(v => d3.interpolateOrRd(logScale(v)));

    const rtAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([colorDomain[0], colorDomain[2]]);

    const caseAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([minCases, maxCases]);

    const pExceedAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain(pExceedColorDomain[0], pExceedColorDomain[2]);

    const rtAxisFn = () => d3.axisBottom(rtAxisScale)
        .tickValues(rtColorScale.ticks())
        .tickFormat(rtColorScale.tickFormat())
        .tickSize(-barHeight);

    const caseAxisFn = () => d3.axisBottom(caseAxisScale)
        .tickValues(logScale.ticks(2))
        .tickFormat(d => d)
        .tickSize(-barHeight);

    const pexceedAxisFn = () => d3.axisBottom(pExceedAxisScale)
        .tickValues(logScale.ticks(2))
        .tickFormat(d => d)
        .tickSize(-barHeight);  

    const axisBottom = map_svg.append("g")
        .attr("class", `x-axis`)
        .attr("transform", `translate(${width - barHeight},${margin.top}) rotate(90)`)
        .call(rtAxisFn)
        .selectAll("text")
        .attr("transform", "translate(-5, 15) rotate(-90)");
    
    const availableDates = rtData.get(rtData.keys()[0]).map(r=>r.Date);
    const bisectDate = d3.bisector(d=>d).left;

    rtFillFn = date => {
        const idx = bisectDate(availableDates, date, 1);
        return d => {  // Fill based on value of Rt
            const rtSeries = rtData.get(d.properties.lad20nm);
            if (!rtSeries) {
                return "#ccc";
            }
            const rt = rtSeries[idx];
            return rtColorScale(rt.Rt50);
        }
    }

    pExceedFillFn = date => {
        const idx = bisectDate(availableDates, date, 1);
        return d => {  // Fill based on value of PExceed
            const pExceedSeries = pexceedData.get(d.properties.lad20nm);
            if (!pExceedSeries) {
                return "#ccc";
            }
            const pExceed = pExceedSeries[idx];
            return pExceedColorScale(pExceed.P10);
        }
    }

    caseFillFn = d => { // Fill based on value of case projection
        const caseProj = nextWeekCaseProjPer100k.get(d.properties.lad20nm);
        if (!caseProj) {
            return "#ccc";
        }
        return caseColorScale(caseProj.caseProjMedian);
    }

    // Draw the map
    const map = g.selectAll("path")
        .data(topojson.feature(topo, topo.objects.Local_Authority_Districts__May_2020__Boundaries_UK_BFC).features)
        .enter().append("path")
        .attr("fill", rtFillFn(d3.max(availableDates)))
        .attr("id", d => "path-" + areaNameToId(d.properties.lad20nm))
        .style("fill-opacity", 1)
        .on("mouseover", function (d) {  // Add Tooltip on hover
            tooltip_div.transition()
                .duration(200)
                .style("opacity", .9);

            tooltip_header.text(d.properties.lad20nm);
            tooltip_info1.text(`Last week Rt: ${getRtForArea(d.properties.lad20nm)}`);
            tooltip_info2.text(`Last week cases (per 100k): ${getCaseHistoryPer100kForArea(d.properties.lad20nm).casesLast7Day}`);
            tooltip_info3.text(`Next week cases (per 100k): ${getCaseProjPer100kForArea(d.properties.lad20nm)}`);
            tooltip_info4.text(`P(Rt > 1): ${getPExceedForArea(d.properties.lad20nm)}`);

            tooltip_div
                .style("left", (d3.event.pageX + 20) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
            d3.select(this).style("fill-opacity", 0.5);
        })
        .on("mouseout", function (d) {
            tooltip_div.style("opacity", 0);
            d3.select(this).style("fill-opacity", 1);
        })
        .on("click", function (d) {
            selectArea(d.properties.lad20nm);
        })
        .attr("d", path)
        .attr("class", "feature")
    map_svg.transition().duration(500).call(zoom.scaleBy, 0.5);

    g.append("path")
        .datum(topojson.mesh(topo, topo.objects.Local_Authority_Districts__May_2020__Boundaries_UK_BFC, (a, b) => a !== b ))
        .attr("class", "mesh")
        .attr("d", path);

    // Draw the color scale
    const defs = map_svg.append("defs");

    const rtGradient = defs.append("linearGradient")
        .attr("id", "rt-gradient");

    const caseGradient = defs.append("linearGradient")
        .attr("id", "case-gradient");

    const pExceedGradient = defs.append("linearGradient")
        .attr("id", "pexceed-gradient");

    rtGradient.selectAll("stop")
        .data([0.5,0.75,.9,1.0,1.1,1.5,2.0].map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: rtColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    caseGradient.selectAll("stop")
        .data(logScale.ticks().map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: caseColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    pExceedGradient.selectAll("stop")
        .data([0.0, 0.25, 0.5, 0.75, 1.0].map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: pExceedColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    const legend = map_svg.append('g')
        .attr("transform", `translate(${width},${barHeight}) rotate(90)`)
        .append("rect")
        .attr('transform', `translate(${margin.left}, 0)`)
        .attr("width", barWidth)
        .attr("height", barHeight)
        .style("fill", "url(#rt-gradient)");

    legendText = map_svg.append("text")
        .attr("x", width - barHeight / 2)
        .attr("y", margin.top + 30)
        .attr("text-anchor", "middle")
        .style("font-size", "12px")
        .text("Rt");

    legendAxis = d3.axisLeft(
        d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([colorDomain[0], colorDomain[2]])
    );
        
    legendBar = map_svg.append("g")
        .attr("transform", `translate(${width - barHeight},${margin.top})`)
        .attr("width", barHeight)
        .attr("height", barWidth)
        .call(legendAxis);

    // Add Rt vs case projection selection
    const showRt = map_svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top + 10)
        .style("font-size", "16px")
        .style("cursor", "pointer")
        .attr("class", "active")
        .text("Rt");

    const showCases = map_svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top + 30)
        .style("font-size", "16px")
        .style("cursor", "pointer")
        .text("Case Projections (Per 100k)");

    const showPExceed = map_svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top + 50)
        .style("font-size", "16px")
        .style("cursor", "pointer")
        .text("P(Rt > 1)");
    
    showRt.on("click", () => {
        if (showCases.classed("active") || showPExceed.classed("active")) {
            map.attr("fill", rtFillFn(d3.max(availableDates)));
            legend.style("fill", "url(#rt-gradient)");
            legendBar.call(d3.axisLeft(
                d3.scaleLinear()
                    .range([margin.left, margin.left + barWidth])
                    .domain([colorDomain[0], colorDomain[2]]))
            );
            axisBottom.call(rtAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("Rt");

            sliderSvg.style("visibility", "visible");
            changeSlider(rtSliderChangeFn);
        }
        showRt.attr("class", "active");
        showCases.attr("class", "");
        showPExceed.attr("class", "");
    });

    showCases.on("click", () => {
        if (showRt.classed("active") || showPExceed.classed("active")) {
            map.attr("fill", caseFillFn);
            legend.style("fill", "url(#case-gradient)");
            legendBar.call(d3.axisLeft(caseLogScale).tickFormat(d3.format(",.0f")));
            axisBottom.call(caseAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("Cases");
            sliderSvg.style("visibility", "hidden");
        }
        showRt.attr("class", "");
        showCases.attr("class", "active");
        showPExceed.attr("class", "");
    });

    showPExceed.on("click", () => {
        if (showRt.classed("active") || showCases.classed("active")) {
            map.attr("fill", pExceedFillFn(d3.max(availableDates)));
            legend.style("fill", "url(#pexceed-gradient)");
            legendBar.call(d3.axisLeft(
                d3.scaleLinear()
                    .range([margin.left, margin.left + barWidth])
                    .domain([pExceedColorDomain[0], pExceedColorDomain[2]]))
            );
            axisBottom.call(pexceedAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("P(Rt > 1)");

            sliderSvg.style("visibility", "visible");
            changeSlider(pExceedSliderChangeFn);
        }
        showCases.attr("class", "");
        showRt.attr("class", "");
        showPExceed.attr("class", "active");
    })

    rtSliderChangeFn = (date) => {
        sliderValueLabel.text(dateFormat(date));
        map.transition().duration(50).attr("fill", rtFillFn(date));
    }

    pExceedSliderChangeFn = (date) => {
        sliderValueLabel.text(dateFormat(date));
        map.transition().duration(50).attr("fill", pExceedFillFn(date));
    }

    const timeSlider = d3.sliderHorizontal()
        .ticks(5)
        .min(d3.min(availableDates))
        .max(d3.max(availableDates))
        .marks(availableDates)
        .tickFormat(d3.timeFormat("%b"))
        .width(sliderWidth)
        .displayValue(false)
        .value(d3.max(availableDates))
        .on('onchange', rtSliderChangeFn);
    
    sliderG.call(timeSlider);
    sliderValueLabel.text(dateFormat(d3.max(availableDates)));

    changeSlider = (changeFn) => {
        timeSlider.min(d3.min(availableDates))
            .on('onchange', changeFn)
            .value(d3.max(availableDates));
            
        sliderValueLabel.text(dateFormat(d3.max(availableDates)));
    }

    // Set up local area search box
    const areas = rtData.keys();

    new autoComplete({
        selector: '#areaSearch',
        minChars: 1,
        source: function(term, suggest){
            term = term.toLowerCase();
            var matches = [];
            for (let i = 0; i < areas.length; i++) {
                if (~areas[i].toLowerCase().indexOf(term)) {
                    matches.push(areas[i]);
                }
            }
            suggest(matches);
        },
        onSelect: (e, term, item) => {
            selectArea(term);
            document.getElementById("#areaSearch").value = "";
        }
    });
}

function plotCaseChart(chartData, projectionData, predictionData, area) {
    const xDomain = d3.extent([...chartData.map(c => c.Date), ...projectionData.map(p => p.Date)]);
    const yDomain = [0, d3.max([...chartData.map(c=>c.cases_new), ...projectionData.map(p=>p.C_median)])];

    const x = d3.scaleTime()
        .domain(xDomain)
        .range([chartMargin.left, chartMargin.left + chartWidth]);
    const y = d3.scaleLinear()
        .domain(yDomain)
        .range([chartHeight, 0]);

    // Define the lines
    const actualCasesLine = d3.line()
        .x(function (d) { return x(d.Date); })
        .y(function (d) { return y(d.cases_new); });

    const smoothedCasesLine = d3.line()
        .x(function (d) { return x(d.Date); })
        .y(function (d) { return y(d.cases_new_smoothed); });

    const predictedCasesLine = d3.line()
        .x(function (d) { return x(d.Date); })
        .y(function (d) { return y(d.C_median); });

    const predictedCasesArea = d3.area()
        .x(function (d) { return x(d.Date); })
        .y0(function (d) { return y(d.C_lower); })
        .y1(function (d) { return y(d.C_upper); });


    const projectedCasesLine = d3.line()
        .x(function (d) { return x(d.Date); })
        .y(function (d) { return y(d.C_median); });

    const projectedCasesArea = d3.area()
        .x(function (d) { return x(d.Date); })
        .y0(function (d) { return y(d.C_lower); })
        .y1(function (d) { return y(d.C_upper); });

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

    predictedArea
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedCasesArea);

    predictedChartLine
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedCasesLine);


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

    caseChartXAxis.call(d3.axisBottom(x));
    caseChartYAxis.call(d3.axisLeft(y).ticks(5));
    caseChartTitle.text(`COVID-19 Cases for ${area}`);

    const predictionDate = d3.min(predictionData.map(c => c.Date));
    const projectionDate = d3.max(chartData.map(c => c.Date));

    const focus = caseChartSvg.append("g")
        .attr("class", "focus")
        .style("display", "none");

    focus.append("line")
        .attr("class", "x-hover-line hover-line")
        .attr("y1", 0)
        .attr("y2", chartHeight);

    focus.append("circle")
        .attr("r", 2);

    focus.append("text")
        .attr("x", 15)
        .attr("dy", ".31em");

    caseChartSvg.append("rect")
        .attr("transform", "translate(" + chartMargin.left + ",0)")
        .attr("class", "overlay")
        .attr("width", chartWidth)
        .attr("height", chartHeight)
        .on("mouseover", function () { focus.style("display", null); })
        .on("mouseout", function () { focus.style("display", "none"); })
        .on("mousemove", mousemove);

    const bisectDate = d3.bisector(function (d) { return d.Date; }).left;
    const allData = [...chartData, ...projectionData];

    function getValue(d) {
        if (d.Date > projectionDate) {
            return d.C_median;
        }
        else {
            return d.cases_new;
        }
    }

    function mousemove() {
        const x0 = x.invert(d3.mouse(this)[0] + chartMargin.left),
            i = bisectDate(allData, x0, 1),
            d0 = allData[i - 1],
            d1 = allData[i],
            d = x0 - d0.Date > d1.Date - x0 ? d1 : d0;

        // TODO: Change this to cases actual, smoothed and projections
        focus.attr("transform", "translate(" + x(d.Date) + "," + y(getValue(d)) + ")");
        focus.select("text").text(function () { return Math.round(getValue(d)); });
        focus.select(".x-hover-line").attr("y2", chartHeight - y(getValue(d)));
    }
}

function plotRtChart(rtData, chartData, projectionData, predictionData, area) {
    const xDomain = d3.extent([...chartData.map(c => c.Date), ...projectionData.map(p => p.Date)]);
    //const xDomain = d3.extent(rtData.map(c => c.Date));
    const yDomain = [0, 3.1];

    const x = d3.scaleTime()
        .domain(xDomain)
        .range([chartMargin.left, chartMargin.left + chartWidth]);
    const y = d3.scaleLinear()
        .domain(yDomain)
        .range([chartHeight, 0]);

    // Define the lines
    const rtMedianLine = d3.line()
        .x(function (d) { return x(d.Date); })
        .y(function (d) { return y(d.Rt50); });

    const rtInnerArea = d3.area()
        .x(function (d) { return x(d.Date); })
        .y0(function (d) { return y(d.Rt25); })
        .y1(function (d) { return y(d.Rt75); });
    
    const rtOuterArea = d3.area()
        .x(function (d) { return x(d.Date); })
        .y0(function (d) { return y(d.Rt025); })
        .y1(function (d) { return y(d.Rt975); });

    rtChartMedianLine
        .datum(rtData)
        .transition()
        .duration(500)
        .attr("d", rtMedianLine);

    rtChartInnerArea
        .datum(rtData)
        .transition()
        .duration(500)
        .attr("d", rtInnerArea);

    rtChartOuterArea
        .datum(rtData)
        .transition()
        .duration(500)
        .attr("d", rtOuterArea);

    rtHorizontalLine
        .attr("transform", `translate(${chartMargin.left}, ${y(1.0)})`)
        .append("line")
        .attr("x2", chartWidth);

    rtChartXAxis.call(d3.axisBottom(x));
    rtChartYAxis.call(d3.axisLeft(y).ticks(5));
    rtChartTitle.text(`Estimated Rt ${area}`);

    //TODO: Refactor this out into a function!
    const focus = rtChartSvg.append("g")
        .attr("class", "focus")
        .style("display", "none");

    focus.append("line")
        .attr("class", "x-hover-line hover-line")
        .attr("y1", 0)
        .attr("y2", chartHeight);

    focus.append("circle")
        .attr("r", 2);

    focus.append("text")
        .attr("x", 15)
        .attr("dy", ".31em");

    rtChartSvg.append("rect")
        .attr("transform", "translate(" + chartMargin.left + ",0)")
        .attr("class", "overlay")
        .attr("width", chartWidth)
        .attr("height", chartHeight)
        .on("mouseover", function () { focus.style("display", null); })
        .on("mouseout", function () { focus.style("display", "none"); })
        .on("mousemove", mousemove);

    const bisectDate = d3.bisector(function (d) { return d.Date; }).left;

    function mousemove() {
        const x0 = x.invert(d3.mouse(this)[0] + chartMargin.left),
            i = bisectDate(rtData, x0, 1),
            d0 = rtData[i - 1],
            d1 = rtData[i],
            d = x0 - d0.Date > d1.Date - x0 ? d1 : d0;

        focus.attr("transform", "translate(" + x(d.Date) + "," + y(d.Rt50) + ")");
        focus.select("text").text(function () { return Math.round(d.Rt50); });
        focus.select(".x-hover-line").attr("y2", chartHeight - y(d.Rt50));
    }
}

function selectArea(selectedArea) {
    let area = selectedArea;
    if (groupedAreaMap.has(selectedArea)) {
        area = groupedAreaMap.get(selectedArea);
        const otherAreas = groupedAreaConstituents.get(area).join(", ");
        d3.select("#sub-heading").text(`Data shown for ${area}, including ${otherAreas}`);
        //d3.select("#cases-title").text(`Cases for ${area} (including ${selectedArea})`);
        //d3.select("#estimates-title").text(`Estimates for ${area} (including ${selectedArea})`);
    }
    else {
        d3.select("#sub-heading").text(`Data shown for ${area}`);
        //d3.select("#cases-title").text(`Cases`);
        //d3.select("#estimates-title").text(`Estimates`);
    }

    g.selectAll("path")
        .style("zIndex", "inherit")
        .style("vector-effect", "non-scaling-stroke")
        .style("stroke-linecap", "round")
        .style("stroke-linejoin", "round")
        .style("stroke-width", "0.5px")
        .style("stroke","#fff");

    d3.select("#path-" + areaNameToId(selectedArea))
        .style("zIndex", "1")
        .style("vector-effect", "non-scaling-stroke")
        .style("stroke-linecap", "round")
        .style("stroke-linejoin", "round")
        .style("stroke-width", "2px")
        .style("stroke", "#222");

    d3.select("#data-heading").text(selectedArea);

    const chartData = caseTimeseries.get(area);
    if (!chartData) {
        console.log("ERROR: No chart data found for area ", area);
        return;
    }
    const projectionData = caseProjTimeseries.get(area);
    if (!projectionData) {
        console.log("ERROR: No projection data found for area ", area);
        return;
    }
    const predictionData = casePredTimeseries.get(area);
    if (!predictionData) {
        console.log("ERROR: No prediction data found for area ", area);
        return;
    }


    const rtChartData = rtData.get(area);
    if (!rtChartData) {
        console.log("ERROR: No Rt data found for area ", area);
        return;
    }

    plotCaseChart(chartData, projectionData, predictionData, area);
    plotRtChart(rtChartData, chartData, projectionData, predictionData, area);

    const caseHistory = getCaseHistoryForArea(area);
    casesLast7Info.text(caseHistory.casesLast7Day);
    casesTotalInfo.text(caseHistory.casesTotal);
    const caseHistoryPer100k = getCaseHistoryPer100kForArea(area);
    casesLast7PerInfo.text(caseHistoryPer100k.casesLast7Day);
    rtInfo.text(getRtForArea(area));
    caseProjInfo.text(getCaseProjForArea(area));
    caseProjPer100kInfo.text(getCaseProjPer100kForArea(area));
}

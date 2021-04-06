// Constants
const TOPOJSON_PATH = "assets/data/uk_lad_boundaries.json";
const SITE_DATA_PATH = "assets/data/site_data.csv";
const NHS_SCOTLAND_MAP = "assets/data/nhs_scotland_health_boards.csv";
const ENGLAND_META_AREA_MAP = "assets/data/england_meta_areas.csv";
const METADATA_PATH = "assets/data/metadata.csv";
const MAP_PATH = 'default';
const RT_PATH = "Rt.csv";
const CASE_PROJECTION_PATH = "Cproj.csv";
const CASE_PREDICTION_PATH = "Cpred.csv";
const INFECTION_PROJECTION_PATH = "Xproj.csv";
const INFECTION_PREDICTION_PATH = "Xpred.csv";
const CASE_WEEKLY_PATH = "Cweekly.csv";
const PEXCEED_PATH = "Pexceed.csv";

const MONTHS = ["Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec"];

const map_svg = d3.select("#map");

const barHeight = 20;
const sliderWidth = 260;
const barWidth = 180;
const margin = ({ top: 20, right: 40, bottom: 30, left: 40 });
const g = map_svg.append("g");

const sliderLeft = 75;
const sliderSvg = d3.select("#slider-svg").style("display", "block");
const sliderG = sliderSvg.append("g")
    .attr('transform', `translate(${sliderLeft},10)`);
sliderSvg.append("text")
    .attr("x", sliderLeft-60)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .style("font-size", "13px")
    .text("Date");
const sliderValueLabel = sliderSvg.append("text")
    .attr("x", sliderLeft + sliderWidth + 85)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .style("font-size", "13px");

const sliderLeftG = sliderSvg.append("g")
    .attr('transform', `translate(${sliderLeft-35},0)`)
    .attr("fill", "none")
    .attr("stroke", "currentColor");

const sliderRightG = sliderSvg.append("g")
    .attr('transform', `translate(${sliderLeft + sliderWidth + 17},0)`)
    .attr("fill", "none")
    .attr("stroke", "currentColor");

const sliderLeftRect = sliderLeftG.append("rect")
    .attr("class", "slider-button");

const sliderLeftPath = sliderLeftG.append("path")
    .attr("transform", "translate(-2, -2)")
    .attr("d", "M15 19l-7-7 7-7")
    .attr("class", "slider-button-path");

const sliderRightRect = sliderRightG.append("rect")
    .attr("class", "slider-button");

const sliderRightPath = sliderRightG.append("path")
    .attr("d", "M9 5l7 7-7 7")
    .attr('transform', `translate(-2, -2)`)
    .attr("class", "slider-button-path");
    

// Set up dimensions for chart
const chartMargin = { top: 30, right: 30, bottom: 30, left: 30 };


const caseChartSvg = d3.select("#case-chart");
// Note: Assuming the two charts are the same size! 
const chartWidth = Math.max(500+caseChartSvg.attr("width") - chartMargin.left - chartMargin.right, 0);
const chartHeight = Math.max(200+caseChartSvg.attr("height") - chartMargin.top - chartMargin.bottom, 0);

const chartCG = caseChartSvg.append("g")
    .attr("transform", "translate(" + chartMargin.left + "," + chartMargin.top + ")");

const actualCChartLine = caseChartSvg.append("path")
    .attr("class", "actual-cases-line")

const smoothedCChartLine = caseChartSvg.append("path")
    .attr("class", "smoothed-cases-line")

const predictedCChartLine = caseChartSvg.append("path")
    .attr("class", "predicted-cases-median-line");

const predictedCInnerArea = caseChartSvg.append("path")
    .attr("class", "predicted-cases-inner-area");

const predictedCOuterArea = caseChartSvg.append("path")
    .attr("class", "predicted-cases-outer-area");

const projectedCChartLine = caseChartSvg.append("path")
    .attr("class", "projected-cases-median-line");

const projectedCInnerArea = caseChartSvg.append("path")
    .attr("class", "projected-cases-inner-area");

const projectedCOuterArea = caseChartSvg.append("path")
    .attr("class", "projected-cases-outer-area");

// Note: Assuming the two charts are the same size! 
const infectionChartSvg = d3.select("#infection-chart");

const chartXG = infectionChartSvg.append("g")
    .attr("transform", "translate(" + chartMargin.left + "," + chartMargin.top + ")");

const actualXChartLine = infectionChartSvg.append("path")
    .attr("class", "actual-cases-line")

const smoothedXChartLine = infectionChartSvg.append("path")
    .attr("class", "smoothed-cases-line")


const predictedXChartLine = infectionChartSvg.append("path")
    .attr("class", "predicted-infections-median-line");

const predictedXInnerArea = infectionChartSvg.append("path")
    .attr("class", "predicted-infections-inner-area");

const predictedXOuterArea = infectionChartSvg.append("path")
    .attr("class", "predicted-infections-outer-area");

const projectedXChartLine = infectionChartSvg.append("path")
    .attr("class", "projected-infections-median-line");

const projectedXInnerArea = infectionChartSvg.append("path")
    .attr("class", "projected-infections-inner-area");

const projectedXOuterArea = infectionChartSvg.append("path")
    .attr("class", "projected-infections-outer-area");


const caseCurrentDateLine = caseChartSvg.append("g");

caseCurrentDateLine.append("line")
    .attr("class", "current-date-line")
    .attr("y1", 0)
    .attr("y2", chartHeight);

const infectionCurrentDateLine = infectionChartSvg.append("g");

infectionCurrentDateLine.append("line")
    .attr("class", "current-date-line")
    .attr("y1", 0)
    .attr("y2", chartHeight);


const rtChartSvg = d3.select("#rt-chart");

const rtChartMedianLine = rtChartSvg.append("path")
    .attr("class", "rt-line");

const rtChartInnerArea = rtChartSvg.append("path")
    .attr("class", "rt-inner-area");

const rtChartOuterArea = rtChartSvg.append("path")
    .attr("class", "rt-outer-area");

const rtProjectedMedianLine = rtChartSvg.append("path")
    .attr("class", "rt-projected-line");

const rtProjectedInnerArea = rtChartSvg.append("path")
    .attr("class", "rt-projected-inner-area");

const rtProjectedOuterArea = rtChartSvg.append("path")
    .attr("class", "rt-projected-outer-area");

const rtHorizontalLine = rtChartSvg.append("g")
    .attr("class", "rt-horizontal-line");

const rtCurrentDateLine = rtChartSvg.append("g");

rtCurrentDateLine.append("line")
    .attr("class", "current-date-line")
    .attr("y1", 0)
    .attr("y2", chartHeight);

// Add the X Axes
const caseChartXAxis = caseChartSvg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chartHeight})`);

const infectionChartXAxis = infectionChartSvg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chartHeight})`);


const rtChartXAxis = rtChartSvg.append("g")
    .attr("class", "chart-x-axis")
    .attr("transform", `translate(0,${chartHeight})`);

const caseX = d3.scaleTime()
    .range([chartMargin.left, chartMargin.left + chartWidth]);

const infectionX = d3.scaleTime()
    .range([chartMargin.left, chartMargin.left + chartWidth]);

const rtX = d3.scaleTime()
    .range([chartMargin.left, chartMargin.left + chartWidth]);

// Add the Y Axis
const caseChartYAxis = caseChartSvg.append("g")
    .attr("class", "chart-y-axis")
    .attr("transform", `translate(${chartMargin.left},0)`);

const infectionChartYAxis = infectionChartSvg.append("g")
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

const infectionChartTitle = infectionChartSvg.append("text")
    .attr("x", (chartWidth / 2))
    .attr("y", (chartMargin.top / 2))
    .attr("text-anchor", "middle")
    .style("font-size", "16px");


const rtChartTitle = rtChartSvg.append("text")
    .attr("x", (chartWidth / 2))
    .attr("y", (chartMargin.top / 2))
    .attr("text-anchor", "middle")
    .style("font-size", "16px");

const casesStartDateInfo = d3.select("#last7-start-date");
const casesEndDateInfo = d3.select("#last7-end-date");
const casesLast7Info = d3.select("#cases-last7-info");
const casesLast7PerInfo = d3.select("#cases-last7-per-info");
const casesTotalInfo = d3.select("#cases-total-info");
const rtInfo = d3.select("#rt-info");
const casesProjStartDateInfo = d3.select("#case-proj-start-date");
const casesProjEndDateInfo = d3.select("#case-proj-end-date");
const caseProjInfo = d3.select("#case-proj-info");
const caseProjPer100kInfo = d3.select("#case-proj-per100k-info");

// Map and projection
const projection = d3.geoMercator()
    .center([-7, 58.4])
    .scale(3000)
    .translate([0, 0]);
const path = d3.geoPath().projection(projection);

// Zooming
const zoomIn = map_svg.append("g").append("text")
    .attr("x", 480)
    .attr("y", 290)
    .attr("width", 20)
    .attr("height", 25)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_in")
    .style("cursor", "pointer")
    .text("+");

const zoomOut = map_svg.append("g").append("text")
    .attr("x", 480)
    .attr("y", 330)
    .attr("width", 20)
    .attr("height", 20)
    .attr("text-anchor", "middle")
    .attr("id", "zoom_out")
    .style("cursor", "pointer")
    .text("\u2013");

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
const tooltipDateFormat = d3.timeFormat("%d %b")
const oneDay = (24*60*60*1000) * 1;
const sixDays = (24*60*60*1000) * 6;
const sevenDays = (24*60*60*1000) * 7;
const eightDays = (24*60*60*1000) * 8;

// Data containers
const rtData = d3.map();
const caseTimeseries = d3.map();          // The real historical cases
const caseProjTimeseries = d3.map();
const casePredTimeseries = d3.map();
const infectionProjTimeseries = d3.map();
const infectionPredTimeseries = d3.map();
const caseWeeklyTimeseries = d3.map();    // Actual and projected cases: plotted.
                                          // TODO: The data in CASE_WEEKLY_PATH give the provenance of the weekly cases as `inferred`
                                          // and `projected`. However, epimap.r takes data from the Count (observed) matrix. The labels
                                          // should be corrected.
const pexceedData = d3.map();
const nextWeekCaseProj = d3.map();
const nextWeekCaseProjPer100k = d3.map();
const caseHistory = d3.map();             // The real historical cases: casesLast7Day is the very last week's cases; the cases not used in modeling.
const caseHistoryPer100k = d3.map();      // The real historical cases
const groupedAreaMap = d3.map();
const groupedAreaConstituents = d3.map();
const populations = d3.map();
let selectedDate = null;

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
tooltip_div.append("br");
const tooltip_info2 = tooltip_div.append("span")
    .attr("class", "info-row2");
tooltip_div.append("br");
const tooltip_info3 = tooltip_div.append("span")
    .attr("class", "info-row3");
tooltip_div.append("br");
const tooltip_info4 = tooltip_div.append("span")
    .attr("class", "info-row4");

// Load external data.
// The last observation date (from site_data.csv) is used to determine when we overlay.
// This date is different from the last prediction date and first projection date.
let lastObservationDate = null;
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
        lastObservationDate = d3.max(cases.slice(-1).map(c => c.Date));
        caseHistory.set(area, {
            casesLast7Day: casesLast7Day,
            casesTotal: casesTotal
        });
    });
});


const urlParams = new URLSearchParams(window.location.search);
const map_path = "assets/data/".concat(urlParams.get("map") || MAP_PATH);
const rt_path = map_path.concat("/", RT_PATH);
const case_projection_path = map_path.concat("/", CASE_PROJECTION_PATH);
const case_prediction_path = map_path.concat("/", CASE_PREDICTION_PATH);
const infection_projection_path = map_path.concat("/", INFECTION_PROJECTION_PATH);
const infection_prediction_path = map_path.concat("/", INFECTION_PREDICTION_PATH);
const case_weekly_path = map_path.concat("/", CASE_WEEKLY_PATH);
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
            Rt975: +d.Rt_97_5,
            provenance: d.provenance
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
            P10: +d.P_10,
            provenance: d.provenance
            // Ignoring other fields, P_08,P_09,P_11,P_12,P_15,P_20 for now
        };
        pexceedData.get(d.area).push(current);
    });
});

// Cproj.csv, starting at the first projection date.
let firstProjectionDate = null;
const loadCaseProjections = d3.csv(case_projection_path).then(data => data.forEach(d => {
    if (!caseProjTimeseries.has(d.area)) {
        caseProjTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.C_025 != undefined) {
      d.C_lower95 = +d.C_025;
      d.C_median = +d.C_50;
      d.C_upper95 = +d.C_975;
    }  

    if (d.C_25 != undefined) {
      d.C_lower50 = +d.C_25;
      d.C_upper50 = +d.C_75; 
    } 

    caseProjTimeseries.get(d.area).push(d);
}))
.then(() => {
    firstProjectionDate = d3.min(caseProjTimeseries.get(caseProjTimeseries.keys()[0]).map(r=>r.Date));
    console.log('First projection date: '+firstProjectionDate);
})
.then(() => {
    caseProjTimeseries.each((projections, area) => {
        let caseProjLower = 0, caseProjMedian = 0, caseProjUpper = 0;
        for (let i = 0; i < 7; i++) {
            caseProjLower += projections[i].C_lower95;
            caseProjMedian += projections[i].C_median;
            caseProjUpper += projections[i].C_upper95;
        }
        nextWeekCaseProj.set(area, {
            caseProjLower: Math.round(caseProjLower),
			// TODO: The data from Cweekly.csv used for the tooltip is rounded down. It means that the
			//       displayed total projected cases could differ by one depending on rouning up or
			//       down. The `floor` ensures that the displayed values matches the tooltip.
            caseProjMedian: Math.floor(caseProjMedian),
            caseProjUpper: Math.round(caseProjUpper)
        });
    });
});

const loadInfectionProjections = d3.csv(infection_projection_path).then(data => data.forEach(d => {
    if (!infectionProjTimeseries.has(d.area)) {
        infectionProjTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.X_025 != undefined) {
      d.X_lower95 = +d.X_025;
      d.X_median = +d.X_50;
      d.X_upper95 = +d.X_975;
    }  

    if (d.X_25 != undefined) {
      d.X_lower50 = +d.X_25;
      d.X_upper50 = +d.X_75; 
    } 

    infectionProjTimeseries.get(d.area).push(d);
}))
.then(() => {
    firstProjectionDate = d3.min(infectionProjTimeseries.get(infectionProjTimeseries.keys()[0]).map(r=>r.Date));
    console.log('First projection date: '+firstProjectionDate);
});


// Cpred.csv, ending at the last prediction date.
// This is the last predicted (not projected!) date. Predicted is from observations,
// projected is into the future.
// It should be true that firstProjectionDate = lastPredictionDate + one day.
let lastPredictionDate = null;
const loadCasePredictions = d3.csv(case_prediction_path).then(data => data.forEach(d => {
    if (!casePredTimeseries.has(d.area)) {
        casePredTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.C_025) {
      d.C_lower95 = +d.C_025;
      d.C_median = +d.C_50;
      d.C_upper95 = +d.C_975;
    }  

    if (d.C_25) {
        d.C_lower50 = +d.C_25;
        d.C_upper50 = +d.C_75; 
    }

    casePredTimeseries.get(d.area).push(d);	
}))
.then(() => {
    lastPredictionDate = d3.max(casePredTimeseries.get(casePredTimeseries.keys()[0]).map(r=>r.Date));
    console.log('Last prediction date: '+lastPredictionDate);
});

const loadInfectionPredictions = d3.csv(infection_prediction_path).then(data => data.forEach(d => {
    if (!infectionPredTimeseries.has(d.area)) {
        infectionPredTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    if (d.X_025) {
      d.X_lower95 = +d.X_025;
      d.X_median = +d.X_50;
      d.X_upper95 = +d.X_975;
    }  

    if (d.X_25) {
        d.X_lower50 = +d.X_25;
        d.X_upper50 = +d.X_75; 
    }

    infectionPredTimeseries.get(d.area).push(d);	
}))
.then(() => {
    lastPredictionDate = d3.max(infectionPredTimeseries.get(infectionPredTimeseries.keys()[0]).map(r=>r.Date));
    console.log('Last prediction date: '+lastPredictionDate);
});


const loadCaseWeekly = d3.csv(case_weekly_path).then(data => data.forEach(d => {
    if (!caseWeeklyTimeseries.has(d.area)) {
        caseWeeklyTimeseries.set(d.area, []);
    }
    d.Date = d3.timeParse("%Y-%m-%d")(d.Date);
    d.C_weekly = +d.C_weekly;
    caseWeeklyTimeseries.get(d.area).push(d);
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

const loadCaseWeeklyAndPopulation = Promise.all([
    loadMetadata,
    loadCaseWeekly
]).then(() => {
    caseWeeklyTimeseries.each((series, area) => {
        series.forEach(c => {
            c.C_weekly_100k = (c.C_weekly / getPopulation(area)) * 100000;
        });
    });
})

const casesAndMeta = Promise.all([
    loadCaseProjections,
    loadCasePredictions,
    loadCases,
    loadCaseWeeklyAndPopulation,
    loadNHSScotland,
    loadEnglandMetaAreas
]).then(() => {
    nextWeekCaseProj.each((caseProj, area) => {
        const pop = getPopulation(area) / 100000;
        nextWeekCaseProjPer100k.set(area, {
            caseProjLower: caseProj.caseProjLower / pop,
            caseProjMedian: caseProj.caseProjMedian / pop,
            caseProjUpper: caseProj.caseProjUpper / pop
        });
    });

    caseHistory.each((cases, area) => {
        const pop = getPopulation(area) / 100000;
        caseHistoryPer100k.set(area, {
            casesLast7Day: cases.casesLast7Day / pop,
            casesTotal: cases.casesTotal / pop
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
const bisectDate = d3.bisector(d=>d).left;

let currentSliderChangeFn = null;

function getRtForArea(area, availableDates, someDate) {
    const rtSeries = rtData.get(area);
    if (!rtSeries) {
        return "? [? - ?]";
    }
    let rt = null;
    if (availableDates) {
        const idx = bisectDate(availableDates, someDate, 1);
        rt = rtSeries[idx];
    }
    else {
        [rt] = rtSeries.slice(-1);
    }

    const median = rt.Rt50.toFixed(1);
    const upper = rt.Rt975.toFixed(1);
    const lower = rt.Rt025.toFixed(1);
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

function getPExceedForArea(area, availableDates) {
    if (!pexceedData.has(area)) {
        return "Unknown";
    }
    const pExceedSeries = pexceedData.get(area);

    let pExceed = null;
    if (availableDates) {
        const idx = bisectDate(availableDates, selectedDate, 1);
        pExceed = pExceedSeries[idx];
    }
    else {
        [pExceed] = pExceedSeries.slice(-1)
    }

    const pRtOver1 = pExceed.P10.toFixed(2);
    return `${pRtOver1}`;
}

function getCaseWeeklyForArea(area, availableDates, someDate) {
    if (!caseWeeklyTimeseries.has(area)) {
        return "Unknown";
    }
    const caseWeeklySeries = caseWeeklyTimeseries.get(area);

    let caseWeekly = null;
    if (availableDates) {
        const idx = bisectDate(availableDates, someDate, 1);
        caseWeekly = caseWeeklySeries[idx];
    }
    else {
        [caseWeekly] = caseWeeklySeries.slice(-1);
    }

    const cases = caseWeekly.C_weekly.toFixed(1);
    const cases100k = caseWeekly.C_weekly_100k.toFixed(1);

    return [`${cases}`, `${cases100k}`];
}

function getCaseProjPer100kForArea(area) {
    if (!nextWeekCaseProj.has(area)) {
        return "Unknown";
    }

    const projection = nextWeekCaseProjPer100k.get(area);

    const cprojmedian = projection.caseProjMedian.toFixed(1);
    const cprojlower = projection.caseProjLower.toFixed(1);
    const cprojupper = projection.caseProjUpper.toFixed(1);

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
    return areaName.replace(" and ", "-").replace(" upon ", "-").replace(",", "").replace(" ", "-").replace("'", "");
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
	maxCases = 800;
    maxColorCases = 800;
		
	const caseColorDomain = [0, 200, maxCases];	
	const caseColorScale = d3.scaleDiverging(t => d3.interpolateViridis(1 - t))
        .domain(caseColorDomain);
	
    const rtAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain(colorDomain);

    const caseAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain([minCases, maxCases]);

    const pExceedAxisScale = d3.scaleLinear()
        .range([margin.left, margin.left + barWidth])
        .domain(pExceedColorDomain);

    const rtAxisFn = () => d3.axisBottom(rtAxisScale)
        .tickValues(rtColorScale.ticks())
        .tickFormat(rtColorScale.tickFormat())
        .tickSize(-barHeight);

    const caseAxisFn = () => d3.axisBottom(caseAxisScale)
        .tickValues(caseColorScale.ticks())
        .tickFormat(caseColorScale.tickFormat())
        .tickSize(-barHeight);

    const pexceedAxisFn = () => d3.axisBottom(pExceedAxisScale)
        .tickValues(pExceedColorScale.ticks(2))
        .tickFormat(d => d)
        .tickSize(-barHeight);  

    const axisBottom = map_svg.append("g")
        .attr("class", `x-axis`)
        .attr("transform", `translate(470, ${margin.top}) rotate(90)`)
        .call(rtAxisFn)
        .selectAll("text")
        .attr("transform", "translate(-5, 15) rotate(-90)");
    
    const availableDates = rtData.get(rtData.keys()[0]).map(r=>r.Date);
    selectDate(lastPredictionDate);
    // selectDate(d3.max(availableDates));

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

    caseFillFn = date => { // Fill based on value of case projection
        const idx = bisectDate(availableDates, date, 1);
        return d => {  // Fill based on value of cases
            const caseSeries = caseWeeklyTimeseries.get(d.properties.lad20nm);
            if (!caseSeries) {
                return "#ccc";
            }
            const cases = caseSeries[idx];
            return caseColorScale(cases.C_weekly_100k);
        }
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
            const startDate = new Date();
            // startDate.setTime(selectedDate.getTime() - sevenDays);
			// Comment: The tooltip_info is a bisect and takes the index of the selected (end) date
			// to look up Rt, weekly cases and probability.
			startDate.setTime(selectedDate.getTime() - sixDays);
            tooltip_info1.text(`${tooltipDateFormat(startDate)} - ${tooltipDateFormat(selectedDate)}`);
            tooltip_info2.text(`Rt: ${getRtForArea(d.properties.lad20nm, availableDates, selectedDate)}`);
            tooltip_info3.text(`Cases/100k: ${getCaseWeeklyForArea(d.properties.lad20nm, availableDates, selectedDate)[1]}`);
            tooltip_info4.text(`P(Rt>1): ${getPExceedForArea(d.properties.lad20nm, availableDates)}`);

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
        .data([0.5, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0].map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: rtColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    caseGradient.selectAll("stop")       
		.data([0, 200, 800].map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: caseColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);		

    pExceedGradient.selectAll("stop")
        .data([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].map((t, i, n) => ({ offset: `${100 * i / n.length}%`, color: pExceedColorScale(t) })))
        .enter().append("stop")
        .attr("offset", d => d.offset)
        .attr("stop-color", d => d.color);

    const legend = map_svg.append('g')
        .attr("transform", `translate(500,${barHeight}) rotate(90)`)
        .append("rect")
        .attr('transform', `translate(${margin.left}, 0)`)
        .attr("width", barWidth)
        .attr("height", barHeight)
        .style("fill", "url(#rt-gradient)");

    legendText = map_svg.append("text")
        .attr("x", 480)
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
        .attr("transform", `translate(470,${margin.top})`)
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
        .text("Cases (per 100k)");

    const showPExceed = map_svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top + 50)
        .style("font-size", "16px")
        .style("cursor", "pointer")
        .text("P(Rt > 1)");
    
    showRt.on("click", () => {
        if (showCases.classed("active") || showPExceed.classed("active")) {
            map.attr("fill", rtFillFn(selectedDate));
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

            changeSlider(rtSliderChangeFn);
        }
        showRt.attr("class", "active");
        showCases.attr("class", "");
        showPExceed.attr("class", "");
    });

    showCases.on("click", () => {
        if (showRt.classed("active") || showPExceed.classed("active")) {
            map.attr("fill", caseFillFn(selectedDate));
            legend.style("fill", "url(#case-gradient)");
            legendBar.call(d3.axisLeft(
                d3.scaleLinear()
                    .range([margin.left, margin.left + barWidth])
                    .domain([caseColorDomain[0], caseColorDomain[2]]))
            );
            axisBottom.call(caseAxisFn)
                .selectAll("text")
                .attr("transform", "translate(-5, 15) rotate(-90)");
            legendText.text("Cases");

            changeSlider(caseSliderChangeFn);
        }
        showRt.attr("class", "");
        showCases.attr("class", "active");
        showPExceed.attr("class", "");
    });

    showPExceed.on("click", () => {
        if (showRt.classed("active") || showCases.classed("active")) {
            map.attr("fill", pExceedFillFn(selectedDate));
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

            changeSlider(pExceedSliderChangeFn);
        }
        showCases.attr("class", "");
        showRt.attr("class", "");
        showPExceed.attr("class", "active");
    })

    rtSliderChangeFn = (date) => {
        selectDate(date);
        map.transition().duration(50).attr("fill", rtFillFn(date));
    }

    currentSliderChangeFn = rtSliderChangeFn;

    pExceedSliderChangeFn = (date) => {
        selectDate(date);
        map.transition().duration(50).attr("fill", pExceedFillFn(date));
    }

    caseSliderChangeFn = (date) => {
        selectDate(date);
        map.transition().duration(50).attr("fill", caseFillFn(date));
    }

    const timeSlider = d3.sliderHorizontal()
        .ticks(5)
        .min(d3.min(availableDates))
        .max(d3.max(availableDates))
        .marks(availableDates)
        .tickFormat(d3.timeFormat("%b"))
        .width(sliderWidth)
        .displayValue(false);
    
    sliderG.call(timeSlider);
    sliderValueLabel.text(dateFormat(d3.max(availableDates))); // YW

    changeSlider = (changeFn) => {
        timeSlider
            .on('onchange', _.debounce(changeFn, 100))
            .value(selectedDate);
            
        currentSliderChangeFn = changeFn;
    }

    changeSlider(rtSliderChangeFn);

    sliderLeftG.on("click", () => {
        let updatedDate = new Date();
        updatedDate.setTime(selectedDate.getTime() - sevenDays);
        if (updatedDate < d3.min(availableDates)) {
            updatedDate = d3.min(availableDates);
        }
        timeSlider.value(updatedDate);
        currentSliderChangeFn(updatedDate);
    });

    sliderRightG.on("click", () => {
        let updatedDate = new Date();
        updatedDate.setTime(selectedDate.getTime() + sevenDays);
        if (updatedDate > d3.max(availableDates)) {
            updatedDate = d3.max(availableDates);
        }
        timeSlider.value(updatedDate);
        currentSliderChangeFn(updatedDate);
    });

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
            document.getElementById("areaSearch").value = "";
        }
    });
}

function plotCaseChart(chartData, projectionData, predictionData, area) {
    const xDomain = d3.extent([...chartData.map(c => c.Date), ...projectionData.map(p => p.Date)]);
    const cases_max = d3.max([...chartData.map(c=>c.cases_new)]);
    const projs_max = d3.max([...projectionData.map(p=>p.C_median)]);
    const yDomain = [0, d3.max([cases_max, d3.min([2.0*cases_max, projs_max])])];

    caseX.domain(xDomain);
    // const x = d3.scaleTime()
    //     .domain(xDomain)
    //     .range([chartMargin.left, chartMargin.left + chartWidth]);
    const y = d3.scaleLinear()
        .domain(yDomain)
        .range([chartHeight, 0]);

    // Define the lines
    const actualCasesLine = d3.line()
        .x(function (d) { return caseX(d.Date); })
        .y(function (d) { return y(d.cases_new); });

    const smoothedCasesLine = d3.line()
        .x(function (d) { return caseX(d.Date); })
        .y(function (d) { return y(d.cases_new_smoothed); });

    const predictedCasesLine = d3.line()
        .x(function (d) { return caseX(d.Date); })
        .y(function (d) { return y(d.C_median); });

    const predictedCasesInnerArea = d3.area()
        .x(function (d) { return caseX(d.Date); })
        .y0(function (d) { return y(d.C_lower50); })
        .y1(function (d) { return y(d.C_upper50); });

    const predictedCasesOuterArea = d3.area()
        .x(function (d) { return caseX(d.Date); })
        .y0(function (d) { return y(d.C_lower95); })
        .y1(function (d) { return y(d.C_upper95); });

    const projectedCasesLine = d3.line()
        .x(function (d) { return caseX(d.Date); })
        .y(function (d) { return y(d.C_median); });

    const projectedCasesInnerArea = d3.area()
        .x(function (d) { return caseX(d.Date); })
        .y0(function (d) { return y(d.C_lower50); })
        .y1(function (d) { return y(d.C_upper50); });
    
    const projectedCasesOuterArea = d3.area()
        .x(function (d) { return caseX(d.Date); })
        .y0(function (d) { return y(d.C_lower95); })
        .y1(function (d) { return y(d.C_upper95); });

    actualCChartLine
        .datum(chartData)
        .transition()
        .duration(500)
        .attr("d", actualCasesLine);

    predictedCInnerArea
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedCasesInnerArea);

    predictedCOuterArea
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedCasesOuterArea);

    predictedCChartLine
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedCasesLine);

    projectedCInnerArea
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedCasesInnerArea);

    projectedCOuterArea
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedCasesOuterArea);

    projectedCChartLine
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedCasesLine);

    caseChartXAxis.call(d3.axisBottom(caseX).tickFormat(d3.timeFormat("%b")));
    caseChartYAxis.call(d3.axisLeft(y).ticks(5));
    caseChartTitle.text(`Covid-19 Cases for ${area}`);


    caseCurrentDateLine
        .style("display", null)
        .attr("transform", "translate(" + caseX(selectedDate) + ", 0)");

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
        .attr("x", -30)
        .attr("y", -15)
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
    const allData = [...chartData, ...projectionData.filter(d => (d.Date > lastObservationDate))];

    function getValue(d) {
        if (d.Date > lastObservationDate) {
            return d.C_median;
        }
        else {
            return d.cases_new;
        }
    }

    function mousemove() {
        const x0 = caseX.invert(d3.mouse(this)[0] + chartMargin.left),
            i = bisectDate(allData, x0, 1),
            d0 = allData[i - 1],
     	    d1 = allData[i];

		if (d1 !== undefined) {
			// On dragging a mouse from out of bounds, the bisection may give a date index past the end
			// of allData.
            const d = x0 - d0.Date > d1.Date - x0 ? d1 : d0;

	        // TODO: Change this to cases actual, smoothed and projections
	        var chartDateFormat = d3.timeFormat("%d %b")
	        focus.attr("transform", "translate(" + caseX(d.Date) + "," + y(getValue(d)) + ")");
	        focus.select("text").text(function () { return chartDateFormat(d.Date)+": "+getValue(d); });
	        focus.select(".x-hover-line").attr("y2", chartHeight - y(getValue(d)));
		}
    }
}


function plotInfectionChart(chartData, projectionData, predictionData, area) {
    const xDomain = d3.extent([...chartData.map(c => c.Date), ...projectionData.map(p => p.Date)]);
    const cases_max = d3.max([...chartData.map(c=>c.cases_new)]);
    const projs_max = d3.max([...projectionData.map(p=>p.X_median)]);
    const yDomain = [0, d3.max([cases_max, d3.min([2.0*cases_max, projs_max])])];

    infectionX.domain(xDomain);
    // const x = d3.scaleTime()
    //     .domain(xDomain)
    //     .range([chartMargin.left, chartMargin.left + chartWidth]);
    const y = d3.scaleLinear()
        .domain(yDomain)
        .range([chartHeight, 0]);

    // Define the lines
    const actualCasesLine = d3.line()
        .x(function (d) { return infectionX(d.Date); })
        .y(function (d) { return y(d.cases_new); });

    const predictedInfectionsLine = d3.line()
        .x(function (d) { return infectionX(d.Date); })
        .y(function (d) { return y(d.X_median); });

    const predictedInfectionsInnerArea = d3.area()
        .x(function (d) { return infectionX(d.Date); })
        .y0(function (d) { return y(d.X_lower50); })
        .y1(function (d) { return y(d.X_upper50); });

    const predictedInfectionsOuterArea = d3.area()
        .x(function (d) { return infectionX(d.Date); })
        .y0(function (d) { return y(d.X_lower95); })
        .y1(function (d) { return y(d.X_upper95); });

    const projectedInfectionsLine = d3.line()
        .x(function (d) { return infectionX(d.Date); })
        .y(function (d) { return y(d.X_median); });

    const projectedInfectionsInnerArea = d3.area()
        .x(function (d) { return infectionX(d.Date); })
        .y0(function (d) { return y(d.X_lower50); })
        .y1(function (d) { return y(d.X_upper50); });
    
    const projectedInfectionsOuterArea = d3.area()
        .x(function (d) { return infectionX(d.Date); })
        .y0(function (d) { return y(d.X_lower95); })
        .y1(function (d) { return y(d.X_upper95); });


    actualXChartLine
        .datum(chartData)
        .transition()
        .duration(500)
        .attr("d", actualCasesLine);

    predictedXInnerArea
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedInfectionsInnerArea);

    predictedXOuterArea
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedInfectionsOuterArea);

    predictedXChartLine
        .datum(predictionData)
        .transition()
        .duration(500)
        .attr("d", predictedInfectionsLine);

    projectedXInnerArea
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedInfectionsInnerArea);

    projectedXOuterArea
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedInfectionsOuterArea);
 
    projectedXChartLine
        .datum(projectionData)
        .transition()
        .duration(500)
        .attr("d", projectedInfectionsLine);

    infectionChartXAxis.call(d3.axisBottom(infectionX).tickFormat(d3.timeFormat("%b")));
    infectionChartYAxis.call(d3.axisLeft(y).ticks(5));
    infectionChartTitle.text(`Covid-19 Infections for ${area}`);

    infectionCurrentDateLine
        .style("display", null)
        .attr("transform", "translate(" + infectionX(selectedDate) + ", 0)");

    const focus = infectionChartSvg.append("g")
        .attr("class", "focus")
        .style("display", "none");

    focus.append("line")
        .attr("class", "x-hover-line hover-line")
        .attr("y1", 0)
        .attr("y2", chartHeight);

    focus.append("circle")
        .attr("r", 2);

    focus.append("text")
        .attr("x", -30)
        .attr("y", -15)
        .attr("dy", ".31em");

    infectionChartSvg.append("rect")
        .attr("transform", "translate(" + chartMargin.left + ",0)")
        .attr("class", "overlay")
        .attr("width", chartWidth)
        .attr("height", chartHeight)
        .on("mouseover", function () { focus.style("display", null); })
        .on("mouseout", function () { focus.style("display", "none"); })
        .on("mousemove", mousemove);

    const bisectDate = d3.bisector(function (d) { return d.Date; }).left;
    const allData = [...chartData, ...projectionData.filter(d => (d.Date > lastObservationDate))];

    function getValue(d) {
        if (d.Date > lastObservationDate) {
            return d.X_median;
        }
        else {
            return d.cases_new;
        }
    }

    function mousemove() {
        const x0 = infectionX.invert(d3.mouse(this)[0] + chartMargin.left),
            i = bisectDate(allData, x0, 1),
            d0 = allData[i - 1],
     	    d1 = allData[i];

		if (d1 !== undefined) {
			// On dragging a mouse from out of bounds, the bisection may give a date index past the end
			// of allData.
            const d = x0 - d0.Date > d1.Date - x0 ? d1 : d0;

	        // TODO: Change this to cases actual, smoothed and projections
	        var chartDateFormat = d3.timeFormat("%d %b")
	        focus.attr("transform", "translate(" + infectionX(d.Date) + "," + y(getValue(d)) + ")");
	        focus.select("text").text(function () { return chartDateFormat(d.Date)+": "+getValue(d); });
	        focus.select(".x-hover-line").attr("y2", chartHeight - y(getValue(d)));
		}
    }
}

function plotRtChart(rtData, chartData, projectionData, predictionData, area) {
    const xDomain = d3.extent([...chartData.map(c => c.Date), ...projectionData.map(p => p.Date)]);
    const yDomain = [0, 3.1];

    rtX.domain(xDomain);

    const rtInferredData = rtData.filter(r=>r.provenance === "inferred");
    const rtProjectedData = rtData.filter(r=>r.provenance === "projected");

    const y = d3.scaleLinear()
        .domain(yDomain)
        .range([chartHeight, 0]);

    // Define the lines
    const rtMedianLine = d3.line()
        .x(function (d) { return rtX(d.Date); })
        .y(function (d) { return y(d.Rt50); });

    const rtInnerArea = d3.area()
        .x(function (d) { return rtX(d.Date); })
        .y0(function (d) { return y(d.Rt25); })
        .y1(function (d) { return y(d.Rt75); });
    
    const rtOuterArea = d3.area()
        .x(function (d) { return rtX(d.Date); })
        .y0(function (d) { return y(d.Rt025); })
        .y1(function (d) { return y(d.Rt975); });

    rtChartMedianLine
        .datum(rtInferredData)
        .transition()
        .duration(500)
        .attr("d", rtMedianLine);

    rtChartInnerArea
        .datum(rtInferredData)
        .transition()
        .duration(500)
        .attr("d", rtInnerArea);

    rtChartOuterArea
        .datum(rtInferredData)
        .transition()
        .duration(500)
        .attr("d", rtOuterArea);

    rtProjectedMedianLine
        .datum(rtProjectedData)
        .transition()
        .duration(500)
        .attr("d", rtMedianLine);

    rtProjectedInnerArea
        .datum(rtProjectedData)
        .transition()
        .duration(500)
        .attr("d", rtInnerArea);

    rtProjectedOuterArea
        .datum(rtProjectedData)
        .transition()
        .duration(500)
        .attr("d", rtOuterArea);

    rtHorizontalLine
        .attr("transform", `translate(${chartMargin.left}, ${y(1.0)})`)
        .append("line")
        .attr("x2", chartWidth);

    rtChartXAxis.call(d3.axisBottom(rtX).tickFormat(d3.timeFormat("%b")));
    rtChartYAxis.call(d3.axisLeft(y).ticks(5));
    rtChartTitle.text(`Estimated Rt ${area}`);

    rtCurrentDateLine
        .style("display", null)
        .attr("transform", "translate(" + rtX(selectedDate) + ", 0)");

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
        .attr("x", -30)
        .attr("y", -15)
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
        const x0 = rtX.invert(d3.mouse(this)[0] + chartMargin.left),
            i = bisectDate(rtData, x0, 1),
            d0 = rtData[i - 1],
            d1 = rtData[i];

		if (d1 !== undefined) {
		    // On dragging a mouse from out of bounds, the bisection may give a date index past the end
		    // of allData.
            const d = x0 - d0.Date > d1.Date - x0 ? d1 : d0;

	        var fr = d3.format(".1f");
	        var chartDateFormat = d3.timeFormat("%d %b")
	        focus.attr("transform", "translate(" + rtX(d.Date) + "," + y(d.Rt50) + ")");
	        focus.select("text").text(function () { return chartDateFormat(d.Date) +": "+ fr(d.Rt50); });
	        focus.select(".x-hover-line").attr("y2", chartHeight - y(d.Rt50));
		}
    }
}

function selectArea(selectedArea) {
    let area = selectedArea;
    if (groupedAreaMap.has(selectedArea)) {
        area = groupedAreaMap.get(selectedArea);
        const otherAreas = groupedAreaConstituents.get(area).join(", ");
        d3.select("#sub-heading").text(`Data shown for ${area}, including ${otherAreas}`);
    }
    else {
        d3.select("#sub-heading").text(`Data shown for ${area}`);
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
    const projectionCasesData = caseProjTimeseries.get(area);
    if (!projectionCasesData) {
        console.log("ERROR: No projection cases data found for area ", area);
        return;
    }
    const predictionCasesData = casePredTimeseries.get(area);
    if (!predictionCasesData) {
        console.log("ERROR: No prediction cases data found for area ", area);
        return;
    }
    const projectionInfectionsData = infectionProjTimeseries.get(area);
    if (!projectionInfectionsData) {
        console.log("ERROR: No projection infections found for area ", area);
        return;
    }
    const predictionInfectionsData = infectionPredTimeseries.get(area);
    if (!predictionInfectionsData) {
        console.log("ERROR: No prediction infections found for area ", area);
        return;
    }


    const rtChartData = rtData.get(area);
    if (!rtChartData) {
        console.log("ERROR: No Rt data found for area ", area);
        return;
    }

    plotCaseChart(chartData, projectionCasesData, predictionCasesData, area);
    plotInfectionChart(chartData, projectionInfectionsData, predictionInfectionsData, area);
    plotRtChart(rtChartData, chartData, projectionCasesData, predictionCasesData, area);

	function formatdate(date) {
		return date.getDate()+' '+MONTHS[date.getMonth()]
	}

	// Put last week's dates in website text.
	// This is: 6 days from last prediction day to last prediction day, or
	//          14 days ago to 8 days ago.
	const casesWeekAgoStart = new Date(lastPredictionDate.getTime() - sixDays);		
	const casesWeekAgoEnd = new Date(lastPredictionDate.getTime());
	casesStartDateInfo.text('('+formatdate(casesWeekAgoStart));
	casesEndDateInfo.text(formatdate(casesWeekAgoEnd)+')');

	// Put this week's dates in website text.
	// This is: The first to the 7th projection day, or
	//          a nowcasted projection onto the week that's just passed, i.e. 7 days ago to yesterday.
	// Recall that we project the week that's just passed, as the case numbers are not reliable.
	const projectionThisWeekStart = new Date(lastPredictionDate.getTime() + oneDay);
	const projectionThisWeekEnd = new Date(lastPredictionDate.getTime() + sevenDays);	
	casesProjStartDateInfo.text('('+formatdate(projectionThisWeekStart));
	casesProjEndDateInfo.text(formatdate(projectionThisWeekEnd)+')');

	// The observed cases, 14 days ago to 8 days ago.
	const availableDates = rtData.get(rtData.keys()[0]).map(r=>r.Date);
	const caseHist = getCaseWeeklyForArea(area, availableDates, casesWeekAgoEnd)
	casesLast7Info.text(caseHist[0]);
	casesLast7PerInfo.text(caseHist[1]);

    const caseHistory = getCaseHistoryForArea(area);
    casesTotalInfo.text(caseHistory.casesTotal);	
	
	// The projected cases, 7 days ago to 1 day ago.
    rtInfo.text(getRtForArea(area, availableDates, projectionThisWeekEnd));	

    // The following code will write the median projections, like they appear in the tooltip.
    // const caseProj = getCaseWeeklyForArea(area, availableDates, projectionThisWeekEnd);
    // caseProjInfo.text(caseProj[0]);
    // caseProjPer100kInfo.text(caseProj[1]);

	// The following code adds error bars...
	caseProjInfo.text(getCaseProjForArea(area));
	caseProjPer100kInfo.text(getCaseProjPer100kForArea(area));
}

function selectDate(date) {
    sliderValueLabel.text(dateFormat(date));
    selectedDate = date;
    rtCurrentDateLine
        .transition().duration(100)
        .attr("transform", "translate(" + rtX(selectedDate) + ", 0)");
    caseCurrentDateLine
        .transition().duration(100)
        .attr("transform", "translate(" + caseX(selectedDate) + ", 0)");
    infectionCurrentDateLine
        .transition().duration(100)
        .attr("transform", "translate(" + caseX(selectedDate) + ", 0)");
}

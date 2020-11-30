---
#
# By default, content added below the "---" mark will appear in the home page
# between the top bar and the list of recent posts.
# To change the home page layout, edit the _layouts/home.html file.
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
#
layout: home
---

<head>
    <!-- Load d3.js -->
    <script src="https://d3js.org/d3.v5.js"></script>
    <script src="https://d3js.org/topojson.v1.min.js"></script>	
    <script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>
    <script src="https://d3js.org/d3-geo-projection.v2.min.js"></script>
    <script src="https://unpkg.com/d3-simple-slider"></script>
    <script src="https://cdn.jsdelivr.net/npm/lodash@4.17.20/lodash.min.js"></script>

    <!--<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto:300,300italic,700,700italic">-->
    <!--<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/normalize/8.0.1/normalize.css">-->

    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/milligram/1.4.0/milligram.css">
    <link rel="stylesheet" href="assets/css/main.css"/>

</head>


<!-- This text is above the map. -->
### **Welcome to the UK Local Covid Map!**

This interactive map visualises the historical 
and predicted future developments of the Covid-19 epidemic
across local authorities in the UK.
You'll see statistical estimates for time-varying regional R numbers from case data.

Note that no one knows the *exact R*, and we have to estimate it from data we do have.
Different models to estimate R may therefore give slightly different estimates. Do also look at the larger regional estimates by the
[MRC Biostatistics Unit](https://www.mrc-bsu.cam.ac.uk/tackling-covid-19/nowcasting-and-forecasting-of-covid-19/) 
and the [CMMID](https://epiforecasts.io/covid/posts/national/united-kingdom/).


<p>
<div class="map-container">
<svg id="map" viewBox="0 0 500 400" preserveAspectRatio="xMidYMid meet" border="1px solid black"> </svg>

<div class="row">
<div class="column" width="100%">
<div class="area-search-container">
 <svg class="search-icon" xmlns="http://www.w3.org/2000/svg"
 fill="none" width="24" height="24" stroke="currentColor">
    <path stroke-linecap="round" stroke-linejoin="round"
    stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
    </svg>
<input id="areaSearch" class="search-input" tabindex="1" width="100%" placeholder="Find Local Authority">
</div>
</div>
</div>

<h1 id="data-heading">Select an area in the map...</h1>
<div id="sub-heading">Data shown for England, Wales and Scotland</div>
<div class="row">
    <div class="column">
        <h3 id="cases-title">Cases</h3>
        <div class="info-row">
            <span class="info-heading">Week ending </span>
            <span id="last7-end-date"></span><span class="info-heading">: </span>
            <span id="cases-last7-info"></span>
        </div>
        <div class="info-row">
            <span class="info-heading">Week ending </span>
            <span id="last7-end-date2"></span><span class="info-heading"> per 100k: </span>
            <span id="cases-last7-per-info"></span>
        </div>
        <h3><span class="info-heading">Total cases: </span><span id="cases-total-info"></span></h3>
    </div>
    <div class="column">
        <h3 id="estimates-title">Projected Cases</h3>
        <div class="info-row">
            <span class="info-heading">Week starting </span>
            <span id="case-proj-start-date"></span><span class="info-heading">: </span>						
            <span id="case-proj-info"></span>
        </div>
        <div class="info-row">
            <span class="info-heading">Week starting</span>
            <span id="case-proj-start-date2"></span><span class="info-heading"> per 100k: </span>						
            <span id="case-proj-per100k-info"></span></div>
        <h3><span class="info-heading">Rt: </span><span id="rt-info"></span></h3>
    </div>
</div>
        
<div id="chart-container">
<svg id="chart" viewBox="0 0 500 200"
preserveAspectRatio="xMidYMid meet" ></svg>
</div>

<div id="chart-container">
<svg id="rt-chart" viewBox="0 0 500 200" 
preserveAspectRatio="xMidYMid meet" ></svg>
</div>

</div>
<script src="assets/js/auto-complete.min.js"></script>
<script src="assets/js/map.js"></script>

</p>

<!-- This text is below the map. -->
The R number roughly measures how fast Covid-19 is spreading in society. In the map, we write "Rt" instead of just "R". The "t" indicates "time". We do this because the number is not constant but can go up or down over time, depending on how fast Covid-19 is spreading at a given time.

You can search for or click on a local authority to see its statistics. The blue parts of the graphs show the number of Covid-19 cases historically over time along with the corresponding estimate for Rt. The red parts of the graphs are predictions made by the model.

Definitions for terms used in the map: 
*   **Case** is an infected individual who has tested positive on the given date, 
under either Pillar 1 or Pillar 2 of the UK's testing strategy.
*   **Rt** denotes the reproduction number at a given point in time: how many secondary cases a single primary case will result in on average. 
**Rt** greater than 1 means that the size of the epidemic is increasing exponentially, and less than 1 means it is shrinking. 
*   **Cases (per 100k)** denotes the weekly number of cases per 100,000 population size. For past weeks, this number comes from historical weekly reported number of cases under Pillars 1+2. For future weeks, this is predicted by the model.
*   **P(Rt>1)** denotes the probability (assigned by the model) that Rt is larger than 1 given the observed case counts.



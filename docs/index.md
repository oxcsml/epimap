---
#
# By default, content added below the "---" mark will appear in the home page
# between the top bar and the list of recent posts.
# To change the home page layout, edit the _layouts/home.html file.
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
#
layout: home
title: Home
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

    
    <link rel="stylesheet" href="assets/css/main.css"/>

</head>

<!-- This text is above the map. -->
### **Welcome to the UK Local Covid Map!**

This interactive map visualises the historical 
and predicted future developments of the Covid-19 epidemic
across local authorities in the UK.

We use a [statistical model]({{ site.baseurl }}{% link data-methods.md %}) to estimate time-varying reproduction numbers (R) from UK 
government data on daily reported cases in each lower tier local authority (LTLA). 
Using these estimates, our model also produces predictions of how the number of cases will
rise/fall in the next few weeks. 
These predictions are made by the model assuming no change in the circumstances affecting
transmission rates and testing in each local authority.

We emphasize that no one knows the *true* reproduction numbers, and we estimate them from data we do have using statistical models. 
It is important to treat a model's estimates with caution and be aware of its 
[limitations]({{ site.baseurl }}{% link limitations.md %}).
Different models to estimate R may give slightly different estimates. 
Do look at estimates and predictions given by other groups at
[Cambridge MRC Biostatistics Unit](https://www.mrc-bsu.cam.ac.uk/tackling-covid-19/nowcasting-and-forecasting-of-covid-19/), [LSHTM CMMID](https://epiforecasts.io/covid/posts/national/united-kingdom/) and [Imperial College London](https://imperialcollegelondon.github.io/covid19local/).

This website and the method behind it were developed by a 
[team]({{ site.baseurl }}{% link people.md %}) in the 
[OxCSML](http://csml.stats.ox.ac.uk/) research group at the 
[University of Oxford](https://www.ox.ac.uk)'s 
[Department of Statistics](http://www.stats.ox.ac.uk), 
along with a number of other collaborators.


<div class="map-container">
<div class="row">
<div class="column">
<svg id="slider-svg" viewBox="0 0 500 50" preserveAspectRatio="xMidYMid meet"></svg>
<svg id="map" viewBox="0 0 500 600" preserveAspectRatio="xMidYMid meet"> </svg>
</div>

<div class="column">

<div class="area-search-container">
 <svg class="search-icon" xmlns="http://www.w3.org/2000/svg"
 fill="none" width="24" height="24" stroke="currentColor">
    <path stroke-linecap="round" stroke-linejoin="round"
    stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
</svg>
<input id="areaSearch" class="search-input" tabindex="1" placeholder="Find Local Authority">
</div>

<h1 id="data-heading">Select an area in the map...</h1>
<div id="sub-heading">Data shown for England, Wales and Scotland</div>
<div class="text-row">
    <div class="text-column">
        <h3 id="cases-title">Last week</h3>
        <div class="info-row">
			<span id="last7-start-date"></span><span class="info-heading"> - </span><span id="last7-end-date"></span>
		</div>							
        <div class="info-row">
			<span class="info-heading">Cases: </span><span id="cases-last7-info"></span>
		</div>
        <div class="info-row">
			<span class="info-heading">Cases per 100k: </span><span id="cases-last7-per-info"></span>
		</div>
        <h3><span class="info-heading">Total cases: </span><span id="cases-total-info"></span></h3>		
    </div>
    <div class="text-column">		
        <h3 id="estimates-title">This week</h3>
        <div class="info-row">
			<span id="case-proj-start-date"></span><span class="info-heading"> - </span><span id="case-proj-end-date"></span>
		</div>												
        <div class="info-row">
			<span class="info-heading">Projected cases: </span><span id="case-proj-info"></span>
		</div>
        <div class="info-row">
			<span class="info-heading">Proj cases/100k: </span><span id="case-proj-per100k-info"></span>
		</div>
        <h3><span class="info-heading">Rt: </span><span id="rt-info"></span></h3>
    </div>
</div>

<div id="chart-container">
<svg id="chart" viewBox="0 0 500 200"
preserveAspectRatio="xMidYMid meet" margin-bottom="-2em"></svg>
</div>

<div id="chart-container">
<svg id="rt-chart" viewBox="0 0 500 200" 
preserveAspectRatio="xMidYMid meet" margin-bottom="-2em"></svg>
</div>

</div>
</div>
        
<script src="assets/js/auto-complete.min.js"></script>
<script src="assets/js/map.js"></script>
</div>
<p></p>
<!-- This text is below the map. -->
The R number roughly measures how fast Covid-19 is spreading in society. 
In the map, we write "Rt" instead of just "R". The "t" indicates "time". 
We do this because the number is not constant but can go up or down over time, 
depending on how fast Covid-19 is spreading at a given time.

You can search for or click on a local authority to see its statistics.
The blue parts of the graphs show the number of reported Covid-19 cases historically along with the corresponding estimates for Rt.
The grey parts of the graphs are future predictions made by the model.
**The grey parts include this week**, as the most recent case numbers are incomplete.
For the cases plot, the thin blue line shows the actual number of daily reported cases while the 
thick lines and the light and dark shaded regions show the median and 95% and 50% credible intervals respectively, as predicted by our model. 
For the Rt plot, the thick lines and the light and dark shaded regions show the median and 
95% and 50% credible intervals respectively, for the inferred weekly Rt.

Definitions for terms used on this website: 
* For predictions, the notation X \[Y,Z\] means that the **median** is X while Y and Z gives the lower and upper endpoints of the 95% **credible interval**. That is the model believes that with 95% chance the interval \[Y,Z\] will contain the corresponding quantity, and there is equal chance that the corresponding quantity will be above and below X.
*   **Case** is an infected individual who has tested positive on the given specimen date, 
under either Pillar 1 or Pillar 2 of the UK's testing strategy.
*   **Rt** denotes the reproduction number at a given point in time: how many secondary cases on average will be infected by a single primary case. 
**Rt** greater than 1 means that the size of the epidemic is increasing exponentially, and less than 1 means it is shrinking exponentially. 
*   **Cases (per 100k)** denotes the weekly number of cases per 100,000 population size. For past weeks, this number comes from historical weekly reported number of cases under Pillars 1+2. For future weeks, this is predicted by the model.
*   **P(Rt>1)** denotes the probability (assigned by the model) that Rt is larger than 1 given the observed case counts.



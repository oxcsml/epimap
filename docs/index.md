---
#
# By default, content added below the "---" mark will appear in the home page
# between the top bar and the list of recent posts.
# To change the home page layout, edit the _layouts/home.html file.
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
#
layout: home
---

<style>
.map-container {
    margin: auto;
    width: 95%;
    height: 100%;
    border: 1px solid black;
}

@media (min-width: 40rem) 
{    .map-frame{
        height: 700px !important;
    }    
}
.map-frame{
    width: 100%;
    height: 1500px;
    border: none;
}
</style>

<!-- This text is above the map. -->
### **Welcome to the UK Local Covid Map!**

This interactive map visualises the historical 
and predicted future developments of the Covid-19 epidemic
across local authorities in the UK.
You'll see statistical estimates for time-varying regional R numbers from case data.

Note that no one knows the *exact R*, and we have to estimate it from data we do have. It is important to treat a model's estimates with caution and be aware of its [limitations]({{ site.baseurl }}{% link limitations.md %}).
Different models to estimate R may give slightly different estimates. Do look at the larger regional estimates by the
[MRC Biostatistics Unit](https://www.mrc-bsu.cam.ac.uk/tackling-covid-19/nowcasting-and-forecasting-of-covid-19/), [CMMID](https://epiforecasts.io/covid/posts/national/united-kingdom/) and [Imperial College London](https://imperialcollegelondon.github.io/covid19local/).

<p>
<div class="map-container">
<iframe class="map-frame" src="{{ '/map.html' | prepend: site.baseurl}}" allow="fullscreen">
</iframe>
</div>
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



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
    position: relative;
    width: 1050px;
    height: 675px;
    overflow: visible;
    border: 1px solid black;
    margin: auto;
}
.map-frame{
    position: relative;
    /*
    This height value is a bit of a hack!
    It is there to let the search box overflow into the post
    Not sure what the correct thing to do in this case is
    */
    height: 1000px;
    width: 1100px;
    overflow: visible;
    margin: auto;
    border: 0;
}
</style>

<!-- This text is above the map. -->
Welcome to the UK Local Covid Map. 
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
<iframe class="map-frame" src="{{ '/map.html' | prepend: site.baseurl}}" allow="fullscreen">
</iframe>
</div>
</p>

<!-- This text is below the map. -->
In the map, we write "Rt" instead of just "R". The "t" indicates "time". We do it because the number is not constant but can go up or down over time.
The map's future projections are made with the assumption that Rt stays fixed in the future.

Definitions for terms in the map: 
*   **Case** is an infected individual who has tested positive on the given date, 
under either Pillar 1 or Pillar 2 of the UK's testing strategy.
*   **Rt** denotes the reproduction number: how many secondary cases a single primary case will result in on average. 
**Rt** greater than 1 implies the size of the epidemic is increasing exponentially, and less than 1 means it is shrinking. 
*   **Cases (Per 100k)** denotes either the historical weekly reported number of cases under Pillars 1+2, normalised by population size, or the predicted number in the future weeks.
*   **P(Rt>1)** denotes the probability that Rt is larger than 1 given the observed case counts.


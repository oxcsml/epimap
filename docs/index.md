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
across local authorities in Great Britain.
We hope it may be of use in informing policy makers, 
health protection teams and the public, 
with regards to the state of the epidemic.


<p>
<div class="map-container">
<iframe class="map-frame" src="{{ 'map.html' | prepend: site.baseurl}}" allow="fullscreen">
</iframe>
</div>
</p>

<!-- This text is below the map. -->
Definitions for terms in the map: 
*   **Case** is an infected individual who has tested positive on the given date, 
under either Pillar 1 or Pillar 2 of the UK's testing strategy.
*   **Rt** denotes the reproduction number: how many secondary cases a single primary case will result in on average. 
**Rt** greater than 1 implies the size of the epidemic is increasing exponentially, and less than 1 means it is shrinking. 
*   **Cases (Per 100k)** denotes either the historical weekly reported number of cases under Pillars 1+2, normalised by population size,
the predicted number in the future weeks.
*   **P(Rt>1)** denotes the probability that Rt is larger than 1 given the observed case counts.


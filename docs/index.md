---
#
# By default, content added below the "---" mark will appear in the home page
# between the top bar and the list of recent posts.
# To change the home page layout, edit the _layouts/home.html file.
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
#
layout: home
---


### Map
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

This text is before the map
This text is before the map
This text is before the map
This text is before the map
This text is before the map
This text is before the map
This text is before the map
This text is before the map
This text is before the map
<p>
<div class="map-container">
<iframe class="map-frame" src="{{ 'map.html' | prepend: site.baseurl}}" allow="fullscreen">
</iframe>
</div>
</p>
*   **Rt** denotes the instantaneous reproduction number: how many secondary cases a single primary case will result in on average. If Rt is larger than 1 then size of the Covid-19 pandemic is increasing exponentially. 
*   **Case Projections (Per 100k)** denotes the predicted number of new infections over the next week in a local authority normalised by population size, as inferred from our model.
*   **P(Rt>1)** denotes the probability that Rt is larger than 1 given the observed case counts, as inferred from our model. 
This text is after the map

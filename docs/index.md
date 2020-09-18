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
This text is after the map
This text is after the map
This text is after the map
This text is after the map
This text is after the map
This text is after the map
This text is after the map
This text is after the map

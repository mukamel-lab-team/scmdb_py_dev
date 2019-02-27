<!DOCTYPE html>
<style>

img {
  vertical-align: middle;
  max-height: 450px;
  display: block;
  margin-left: auto;
  margin-right: auto;
}

/* Position the image container (needed to position the left and right arrows) */
.container {
  position: relative;
}

/* Hide the images by default */
.mySlides {
  display: none;
}

/* Add a pointer when hovering over the thumbnail images */
.cursor {
  cursor: pointer;
}

/* Next & previous buttons */
.prev,
.next {
  cursor: pointer;
  position: absolute;
  top: 40%;
  width: auto;
  padding: 16px;
  margin-top: -50px;
  color: black;
  font-weight: bold;
  font-size: 20px;
  border-radius: 0 3px 3px 0;
  user-select: none;
  -webkit-user-select: none;
}

/* Position the "next button" to the right */
.next {
  right: 0;
  border-radius: 3px 0 0 3px;
}

/* On hover, add a black background color with a little bit see-through */
.prev:hover,
.next:hover {
  background-color: rgba(200, 0, 0, 0.8);
}

/* Number text (1/3 etc) */
.numbertext {
  /*color: #f2f2f2;*/
  font-size: 18px;
  padding: 8px 12px;
  position: absolute;
  top: 0;
}

.row:after {
  content: "";
  display: table;
  clear: both;
}

/* Six columns side by side */
.column {
  /*float: left;*/
  width: 16.66%;
  max-height: 100px;
}

/* Add a transparency effect for thumnbail images */
.demo {
  opacity: 0.6;
}

.active,
.demo:hover {
  opacity: 1;
}
</style>

<h2 style="text-align:center">CEMBA - Mouse brain dissections</h2>

<div class="container">

  <?php
    $pattern    = '/var/www/html/CEMBA_dissection_images/CEMBA_Slice*[0-9].png';
    $files = glob($pattern);
    sort($files,SORT_NATURAL);
    foreach ($files as $key => $value) {
      $currfile = basename($value);
      printf('<div class="mySlides">');
      printf(' <div class="numbertext">Slice %d / %d</div>', $key+1, sizeof($files));
      printf(' <img src="https://brainome.ucsd.edu/CEMBA_dissection_images/%s" align="center">',$currfile);
      printf('</div>');
    }
  ?>
  <a class="prev" onclick="plusSlides(-1)">❮ Anterior</a>
  <a class="next" onclick="plusSlides(1)">Posterior ❯</a>

  <div class="row" style="display: inline-block; overflow-x: scroll; overflow-y: hidden; white-space: nowrap; ">
  <?php
    $pattern    = '/var/www/html/CEMBA_dissection_images/CEMBA_Slice*_sm.png';
    $files = glob($pattern);
    sort($files,SORT_NATURAL);
    foreach ($files as $key => $value) {
      $currfile = basename($value);
      printf('<div class="column" style="display: inline-block; overflow-x: scroll; overflow-y: hidden; white-space: nowrap; max-height: 100px;">');
      printf(' <img class="demo cursor" src="https://brainome.ucsd.edu/CEMBA_dissection_images/%s" style="width:100%%" onclick="currentSlide(%d)" alt="Slice %d / %d">', $currfile, $key+1, $key+1, sizeof($files));
      printf('</div>');
    }
  ?>
  </div>
</div>

<script>
var variables;
var queryPost;
var slideIndex = 1;

// Get the slide index to be shown from the URL
function getQueryVariable(variable) 
{
  if (!variables)
  { 
    if (!queryPost) queryPost = window.location.search.substring(1);
    variables = decodeURI(queryPost).split("&");
  }
  for (var i=0;i<variables.length;i++) 
  { 
    var pair = variables[i].split("="); 
    if (pair[0] == variable) return pair[1]; 
  } 
} 
a=getQueryVariable('slideIndex'); if (a) { slideIndex = a;}

showSlides(slideIndex);

function plusSlides(n) {
  showSlides(slideIndex += n);
}

function currentSlide(n) {
  showSlides(slideIndex = n);
}

function showSlides(n) {
  var i;
  var slides = document.getElementsByClassName("mySlides");
  var dots = document.getElementsByClassName("demo");
  // var captionText = document.getElementById("caption");
  if (n > slides.length) {slideIndex = 1}
  if (n < 1) {slideIndex = slides.length}
  for (i = 0; i < slides.length; i++) {
      slides[i].style.display = "none";
  }
  for (i = 0; i < dots.length; i++) {
      dots[i].className = dots[i].className.replace(" active", "");
  }
  slides[slideIndex-1].style.display = "block";
  dots[slideIndex-1].className += " active";
  // captionText.innerHTML = dots[slideIndex-1].alt;
}

currentSlide(slideIndex);

</script>

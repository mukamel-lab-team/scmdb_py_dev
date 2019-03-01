<style>

.myimg {
  vertical-align: middle;
  width: 100%;
  display: block;
  margin-left: auto;
  margin-right: auto;
}

/* Position the image container (needed to position the left and right arrows) */
.container {
  position: relative;
  width: 90%;
  margin: auto;
  overflow: hidden;
}

/* Hide the images by default */
.mySlides {
  display: none;
  margin: left;
  width: 100%;
}

/* Add a pointer when hovering over the thumbnail images */
.cursor {
  cursor: pointer;
}

/* Next & previous buttons */
.prevButton,
.nextButton {
  cursor: pointer;
  position: absolute;
  top: 0px;
  width: auto;
  padding: 16px;
  color: black;
  font-weight: bold;
  font-size: 20px;
  border-radius: 0 3px 3px 0;
  user-select: none;
  -webkit-user-select: none;
}

/* Position the "nextButton button" to the right */
.nextButton {
  right: 0;
  border-radius: 3px 0 0 3px;
}

/* On hover, add a black background color with a little bit see-through */
.prevButton:hover,
.nextButton:hover {
  background-color: #05274A;
  color: white;
}

/* Number text (1/3 etc) */
.numbertext {
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
  width: 33%;
  max-width: 200px;
  display: inline-block;
  white-space: nowrap;
}

.stretch {
  width: 100%;
  display: inline-block;
  font-size: 0;
  line-height: 0
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

<!-- <h2 style="text-align:center">CEMBA - Mouse brain dissections</h2>
 -->
<div class="container">
  <table>
    <tr>
      <td>
        <div style="height:450px; overflow:auto; text-align: justify;">
          <?php
            $pattern    = '/var/www/html/CEMBA_dissection_images/CEMBA_Slice*_sm.png';
            $files = glob($pattern);
            sort($files,SORT_NATURAL);
            foreach ($files as $key => $value) {
              $currfile = basename($value);
              printf('<div class="column">');
              printf(' <img class="demo cursor" src="https://brainome.ucsd.edu/CEMBA_dissection_images/%s" style="width:100%%;" onclick="currentSlide(%d)" alt="Slice %d / %d">', $currfile, $key+1, $key+1, sizeof($files));
              printf('</div>');
            }
          ?>
            <span class="stretch"></span>
          </div>
      </td>
      <td  style="min-width: 600px;">
        <div style="height:450px; overflow:auto;">
          <?php
            $pattern    = '/var/www/html/CEMBA_dissection_images/CEMBA_Slice*[0-9].png';
            $files = glob($pattern);
            sort($files,SORT_NATURAL);
            foreach ($files as $key => $value) {
              $currfile = basename($value);
              $url = sprintf('https://brainome.ucsd.edu/CEMBA_dissection_images/%s', $currfile);
              printf('<div class="mySlides"">');
              printf(' <h2 align="center" style="font-family:inherit;">Slice %d / %d</h2>', $key+1, sizeof($files));
              printf(' <a target="_blank" onclick="updateTable(%d)"><img class="myimg" src="%s" align="center" style="height:400px; width:auto;"></a>', $key+1, $url, $url);
              printf('</div>');
            }
          ?>
          <a class="prevButton" onclick="plusSlides(-1)">❮ Anterior</a>
          <a class="nextButton" onclick="plusSlides(1)">Posterior ❯</a>
        </div>
      </td>
    </tr>
  </table>
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

function updateTable(num) {
  // Fill in the slice number that the user clicked into the search field of the table
  document.getElementsByTagName("input")[0].value = num;
}

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

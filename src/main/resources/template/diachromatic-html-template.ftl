<!doctype html>
<html class="no-js" lang="">

<head>
  <meta charset="utf-8">
  <meta http-equiv="x-ua-compatible" content="ie=edge">
  <title>Diachromatic Analysis</title>
  <meta name="description" content="">
  <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

  <style>
@import url("https://www.jax.org/_res/css/modules/jax-base/p01-fonts.css");
@import url("https://www.jax.org/_res/css/modules/fonts-extended.css");


* {
    -moz-box-sizing: border-box;
    -webkit-box-sizing: border-box;
    box-sizing: border-box
}

html, body, h1, li, a, article, aside, footer, header, main, nav, section {
	padding: 0;
	margin: 0;
}

html, body {
	font-size:14px;
}

body {
	font-family:"DIN Next", Helvetica, Arial, sans-serif;
	line-height:1.25;
	background-color:#e0e3ea;
}


body > header, nav, main, body > section, footer {
max-width:1200px;
margin-left:auto;
margin-right:auto;
}

@media(min-width:1440px) {
body > header, nav, main, body > section, footer {
    width:83.3333%;
    max-width:unset;
    }
}

main, body > section {
	margin-top:1.5rem;
	margin-bottom:1.5rem;
}

body > header, body > section {
	padding:2.1rem 2rem 1.6rem;
}

a[href] {
	color:#05396b;
}

a[href]:hover {
	color:#009ed0;
}

p {
	padding:0;
	margin:0.75rem 0;
}

h1 {
	font-family:"DIN Next", Helvetica, Arial, sans-serif;
	font-weight:700;
	font-size:1.8rem;
	line-height:1;
}

/* Your really should address semantic issues with your markup that make selectors like this necessary */

main > section > a[name="othergenes"] > h3,
h2 {
	font-family:"DIN Next", Helvetica, Arial, sans-serif;
	font-weight:700;
	font-size:1.5rem;
	line-height:1;
	margin:0 0 0.5rem;
	padding:0;
}

h3 {
	font-family:"DIN Next", Helvetica, Arial, sans-serif;
	font-weight:700;
	font-size:1.2rem;
	line-height:1;
	margin:0 0 0.5rem;
	padding:0;
}



main ul, main ol {
	margin:0.5rem 0 0.5rem 1.4rem;
	padding:0;
}

main li {
	margin:0.25rem 0;
	padding:0;
}




.banner {
	background-color: #05396b;
	color: white;
}

nav {
	background-color: #05396b;
	margin-top:1px;
	overflow:auto;
	zoom:1;
	padding:0;
}

nav a[href] {
	color:white;
	text-decoration:none;
	color:rgba(255,255,255,0.8);
	font-size:1.2rem;
	display:block;
	padding:1rem;
	font-weight:400;
}

nav li:last-child a[href] {
	padding-right:2.25rem;
}

nav a[href]:hover {
	color:#05396b;
	background-color:#04c3ff;
}

#navi ul {
	display:table;
	float:right;
	margin:0;
}

#navi li {
	display:block;
	float:left;
}


main > section:first-child {
	margin-top:1.5rem;
	margin-bottom:1.5rem;
	background-color:white;
	padding:2.1rem 2rem 1.6rem;

}

main > section:nth-child(2) {
	margin-top:1.5rem;
	margin-bottom:0;
	background-color:white;
	padding:2.1rem 2rem 1.6rem;

}

main > section + section ~ section > article {
	padding:2.1rem 2rem 1.6rem;
	margin-top:1px;
	background-color:white;
}

main > section > a[name="othergenes"] {
	display:block;
	margin-top:1.5rem;
	background-color:white;
	padding:2.1rem 2rem 1.6rem;
}

table {
	border-collapse: collapse;
	width:100%;
	margin:0.5rem 0;
}

th, td {
	text-align:left;
	padding:0.4rem 0.5rem 0.25rem;
}

th {
	background-color: #e0e3ea;
	border-bottom:1px solid white;
}

table.redTable {
	width:auto;
	min-width:50%;
}

table.redTable td {
	background-color:#f0f3fa;
}

table.minimalistBlack th,
table.minimalistBlack td {
	border:2px solid #e0e3ea;
}

table.minimalistBlack.red td {
	background: red;
}

td.red {
	background-color:#f0f3fa;
}


a[name="othergenes"] table.redTable {

}

a[name="othergenes"] table.redTable td.disease {
	font-size:0.928rem;
	padding-top:0.35rem;
	padding-bottom:0.15rem;
	text-transform: lowercase
}

a[name="othergenes"] table.redTable > tbody > tr:nth-child(even) > td {
	background-color:white;
}

a[name="othergenes"] table.redTable > tbody > tr:hover > td {
	background-color:#cceaff;
}

a[name="othergenes"] table.redTable a {
	text-decoration: none;
	display:block;
}

a[name="othergenes"] table.redTable a:hover {
	text-decoration: underline;
}

a[name="othergenes"] table.redTable a::first-letter {
	text-transform: uppercase;
}

/* Create two equal columns that floats next to each other */
.column {
  float: left;
  width: 50%;
  padding: 10px;
}

/* Clear floats after the columns */
.row:after {
  content: "";
  display: table;
  clear: both;
}


footer {
	background-color: #05396b;
	color: white;
	padding: 1rem 2rem;
}
  </style>
</head>

<body>
  <!--[if lte IE 9]>
    <p class="browserupgrade">You are using an <strong>outdated</strong> browser. Please <a href="https://browsehappy.com/">upgrade your browser</a> to improve your experience and security.</p>
  <![endif]-->
<header class="banner">

<script>${highcharts!""}</script>
<script>${exporting!""}</script>
<script>${exportdata!""}</script>


    <h1>Diachromatic: Preparation and Quality Control of Hi-C and Capture Hi-C Data</h1>
</header>

<nav>
    <div id="navi">
        <ul>
            <li><a href="#sample">Sample</a></li>
            <li><a href="#truncation">Truncation</a></li>
            <li><a href="#align">Alignment</a></li>
             <li><a href="#count">Count</a></li>
            <li><a href="#settings">Settings</a></li>
            <li><a href="#about">About</a></li>
        </ul>
    </div>
</nav>
<main>
  <section>
    <a name="sample"></a>
    <h2>Sample name: </h2>
    <article>
      <div class="row">
        <div class="column" style="background-color:#aaa;">
          <h2>To do show names of files etc.</h2>


    </article>
  </section>


  <section>
   <article>
    <a name="truncation"></a>
    <h2>Truncation statistics</h2>
    <p>
    Truncation was performed with a length threshold of ${length_threshold!"n/a"} nucleotides (reads whose length is below this threshold
    after truncation are removed, denoted <i>Too short to map</i> in the table). A total of ${removed_pairs_one_or_two_reads_too_short!"n/a"}
    pairs contained either a forward or a reverse read (or both) that were too short to map.</p>
    <p>
     <table class="redTable">
         <tr><th></th><th>Forward Read</th><th>Reverse Read</th></tr>
             <tr><td><b>Total Reads</b></td><td>${total_read_pairs_processed!"n/a"}</td><td>${total_read_pairs_processed!"n/a"}</td></tr>
          <tr><td><b>Truncated Reads</b></td><td>${truncated_forward_reads!"n/a"}</td><td>${truncated_reverse_reads!"n/a"}</td></tr>
            <tr><td><b>Dangling Reads</b></td><td>${dangling_forward_reads!"n/a"}</td><td>${dangling_reverse_reads!"n/a"}</td></tr>
     <tr><td><b>Too short to map</b></td><td>${short_removed_forward_reads!"n/a"}</td><td>${short_removed_reverse_reads!"n/a"}</td></tr>
         </table>
    </p>
     <div id="container" style="min-width: 310px; height: 400px; margin: 0 auto"></div>
            <script type="text/javascript">
                ${truncationjs!""}
            		</script>
    </article>
    <#if align??>
    </section>
      <section>
       <article>
        <a name="align"></a>
        <h2>Align statistics</h2>
        <p>
        <ol>
         <#list align as line>
            <li>${line}</li>
         </#list>
        </ol>
        </p>
        </article>
        </section>
    </section>

    </#if>
    <#if count??>
      <section>
       <article>
        <a name="count"></a>
        <h2>Count statistics</h2>
        <p>
        <ol>
         <#list count as line>
            <li>${line}</li>
         </#list>
        </ol>
        </p>

        </article>
        </section>
       </#if> <!-- count section -->
</main>

  <section>
   <article>
    <a name="about"></a>
    <h2>About</h2>

    <p>Diachromatic is a tool for ....</p>
    </article>
    </section>
<footer>
  <p>Diachromatic &copy; 2019</p>
</footer>
</body>
</html>
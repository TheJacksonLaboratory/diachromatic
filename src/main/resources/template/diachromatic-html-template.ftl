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

* {
  box-sizing: border-box;
}

.column {
  float: left;
  width: 50%;
}

/* Clear floats after the columns */
.row:after {
  content: "";
  display: table;
  clear: both;
}

hr {
    display: block;
    height: 1px;
    border: 0;
    border-top: 1px solid #ccc;
    margin: 1em 0;
    padding: 0;
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


<!-- Section on truncation -->
  <section>
    <article>
     <a name="truncation"></a>
     <h2>Truncation statistics</h2>
     <p>Truncation was performed with a length threshold of ${length_threshold!"n/a"} nucleotides (reads whose length is below this threshold
        after truncation are removed, denoted <i>Too short to map</i> in the table). A total of ${removed_pairs_one_or_two_reads_too_short!"n/a"}
         pairs contained either a forward or a reverse read (or both) that were too short to map.
     </p>
     <table class="redTable">
       <tr><th></th><th>Forward Read</th><th>Reverse Read</th></tr>
       <tr><td><b>Total Reads</b></td><td>${total_read_pairs_processed!"n/a"}</td><td>${total_read_pairs_processed!"n/a"}</td></tr>
       <tr><td><b>Truncated Reads</b></td><td>${truncated_forward_reads!"n/a"}</td><td>${truncated_reverse_reads!"n/a"}</td></tr>
       <tr><td><b>Dangling Reads</b></td><td>${dangling_forward_reads!"n/a"}</td><td>${dangling_reverse_reads!"n/a"}</td></tr>
       <tr><td><b>Too short to map</b></td><td>${short_removed_forward_reads!"n/a"}</td><td>${short_removed_reverse_reads!"n/a"}</td></tr>
     </table>
      <div id="container" style="min-width: 310px; height: 400px; margin: 0 auto"></div>
      <script type="text/javascript">${truncationjs!""}</script>
    </article>
  </section>

 <!-- Section on alignment -->

 <section>
   <article>
     <a name="align"></a>
     <h2>Align statistics</h2>
     <p>Diachromatic performs alignment of the truncated reads
        with <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="__blank">bowtie2</a>. TODO-more text.
    </p>
    <div class="row">
      <div class="column" > <!-- upper left quadrant -->
        <h1>Read pair analysis</h1>
        <p>TODO Some text..</p>
        <table class="redTable">
          <tr><th></th><th>Count</th></tr>
          <tr><td><b>Total Reads</b></td><td>${align_total_read_pairs_processed!"n/a"}</td></tr>
          <tr><td><b>Unmapped Read pairs</b></td><td>${align_unmapped_read_pairs!"n/a"}</td></tr>
          <tr><td><b>Multimapped Read Pairs</b></td><td>${align_multimapped_read_pairs!"n/a"}</td></tr>
          <tr><td><b>Paired Read pairs</b></td><td>${align_paired_read_pairs!"n/a"}</td></tr>
          <tr><td><b>Unique paired Read pairs</b></td><td>${align_unique_paired_read_pairs!"n/a"}</td></tr>
          <tr><td><b>Duplicated Read pairss</b></td><td>${align_duplicated_pairs!"n/a"}</td></tr>
         </table>
      </div>
    <div class="column" > <!-- lower right quadrant -->
      <h1>Read  analysis</h1>
      <p>TODO Some text about the reads..</p>
      <table class="redTable">
        <tr><th></th><th>Forward Read</th><th>Reverse Read</th></tr>
        <tr><td><b>Unmapped Reads</b></td><td>${align_unmapped_R1_reads!"n/a"}</td><td>${align_unmapped_R2_reads!"n/a"}</td></tr>
        <tr><td><b>Multimapped Reads</b></td><td>${align_multimapped_R1_reads!"n/a"}</td><td>${align_multimapped_R2_reads!"n/a"}</td></tr>
      </table>
    </div>
   </div>

   <hr/> <!-- dividing line between top two and bottom two quadrants --maybe remove? -->

   <div class="row">
     <div class="column" > <!-- lower left quadrant -->
       <div id="container_alignRead" style="min-width: 150px; height: 200px; margin: 0 auto"></div>
       <script type="text/javascript">${alignjsR!""}</script>
     </div>
     <div class="column" > <!-- lower right quadrant -->
       <div id="container_alignReadPair" style="min-width: 210px; height: 200px; margin: 0 auto"></div>
       <script type="text/javascript">${alignjsRP!""}</script>
     </div>
   </div>

   <p> TODO -- consider reformating the table because many of the items are for read-pairs and not reads</p
   </article>
  </section>



     <section>
           <article>
            <a name="artefact"></a>
            <h2>Artefact statistics</h2>
            <p>
            Artefacts are TODO--Write a summary here
            </p>
            <p>
           <table class="redTable">
             <tr><th>Artefact type</th><th>Count</th></tr>
               <tr><td><b>Unligated</b></td><td>${align_unligated!"n/a"}</td></tr>
               <tr><td><b>unligated_by_size</b></td><td>${align_unligated_by_size!"n/a"}</td></tr>
                <tr><td><b>unligated (same_internal)</b></td><td>${align_unligated_same_internal!"n/a"}</td></tr>
                <tr><td><b>self-ligated</b></td><td>${align_self_ligated!"n/a"}</td></tr>
                <tr><td><b>self-ligated_by_size</b></td><td>${align_self_ligated_by_size!"n/a"}</td></tr>
                  <tr><td><b>self-ligated (same internal)</b></td><td>${align_self_ligated_same_internal!"n/a"}</td></tr>
                   <tr><td><b>chimeric</b></td><td>${align_chimeric!"n/a"}</td></tr>
                    <tr><td><b> chimeric short</b></td><td>${align_chimeric_short!"n/a"}</td></tr>
                     <tr><td><b>chimeric long</b></td><td>${align_chimeric_long!"n/a"}</td></tr>
                      <tr><td><b>chimeric valid</b></td><td>${align_chimeric_valid!"n/a"}</td></tr>
                       <tr><td><b> strange internal </b></td><td>${align_strange_internal!"n/a"}</td></tr>
                        <tr><td><b> dangling_end_pairs_total  </b></td><td>${align_dangling_end_pairs_total!"n/a"}</td></tr>
                         <tr><td><b> trans_pairs_total  </b></td><td>${align_trans_pairs_total!"n/a"}</td></tr>
                </table>
                </p>

                  <p> Detailed results and sanity checks TODO -- better text. </p>
                           <table class="redTable">
                           <caption>dangling end analysis</caption>
                             <tr><th>Artefact type</th><th>Count</th></tr>
                                         <tr><td><b> n_paired_unique_un_ligated_dangling  </b></td><td>${align_n_paired_unique_un_ligated_dangling!"n/a"}</td></tr>
                                         <tr><td><b>  n_paired_unique_self_ligated_dangling </b></td><td>${align_n_paired_unique_self_ligated_dangling!"n/a"}</td></tr>
                                         <tr><td><b> n_paired_unique_too_short_dangling  </b></td><td>${align_n_paired_unique_too_short_dangling!"n/a"}</td></tr>
                                         <tr><td><b> n_paired_unique_too_long_dangling  </b></td><td>${align_n_paired_unique_too_long_dangling!"n/a"}</td></tr>
                                         <tr><td><b> n_paired_unique_valid_dangling  </b></td><td>${align_n_paired_unique_valid_dangling!"n/a"}</td></tr>
                                         <tr><td><b> n_paired_strange_internal_dangling  </b></td><td>${align_n_paired_strange_internal_dangling!"n/a"}</td></tr>
                                </table>
                                <br/> <br/>

                                 <table class="redTable">
                                   <caption>trans pair  analysis</caption>
                                     <tr><th>Artefact type</th><th>Count</th></tr>
                                        <tr><td><b> n_paired_unique_un_ligated_trans  </b></td><td>${align_n_paired_unique_un_ligated_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_paired_unique_self_ligated_trans  </b></td><td>${align_n_paired_unique_self_ligated_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_paired_unique_too_short_trans  </b></td><td>${align_n_paired_unique_too_short_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_paired_unique_too_long_trans  </b></td><td>${align_n_paired_unique_too_long_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_paired_unique_valid_trans  </b></td><td>${align_n_paired_unique_valid_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_paired_strange_internal_trans  </b></td><td>${align_n_paired_strange_internal_trans!"n/a"}</td></tr>
                                        <tr><td><b> n_total_trans  </b></td><td>${align_n_total_trans!"n/a"}</td></tr>
                                 </table>

 <br/> <br/>

                                 <table class="redTable">
                                   <caption>Summary</caption>
                                     <tr><th>Artefact type</th><th>Count</th></tr>
                                        <tr><td><b> Yield of valid pairs (YVP)</b></td><td>${align_YVP!"n/a"}</td></tr>
                                        <tr><td><b> Cross-ligation coefficient (CLC)  </b></td><td>${align_CLC!"n/a"}</td></tr>
                                         <tr><td><b>Re-ligation coefficient (RLC)   </b></td><td>${align_RLC!"n/a"}</td></tr>
                                          <tr><td><b>Hi-C pair duplication rate (HPDR)</b></td><td>${align_HPDR!"n/a"}</td></tr>
                                 </table>

            </article>
            </section>
        </section>



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
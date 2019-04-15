<!doctype html>
<html class="no-js" lang="">
    <head>
        <meta charset="utf-8">
        <meta http-equiv="x-ua-compatible" content="ie=edge">
        <title>Diachromatic Analysis</title>
        <meta name="description" content="">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.1/css/solid.css" integrity="sha384-QokYePQSOwpBDuhlHOsX0ymF6R/vLk/UQVz3WHa6wygxI5oGTmDTv8wahFOSspdm" crossorigin="anonymous">
        <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.1/css/fontawesome.css" integrity="sha384-vd1e11sR28tEK9YANUtpIOdjGW14pS87bUBuOIoBILVWLFnS+MCX9T6MMf0VdPGq" crossorigin="anonymous">
        <style>
            <#include "bootstrap.min.css">
        </style>
        <style>
            <#include "main.css">
        </style>
    </head>
    <body>
    <!--[if lte IE 9]>
    <p class="browserupgrade">You are using an <strong>outdated</strong> browser.
    Please <a href="https://browsehappy.com/">upgrade your browser</a> to improve your experience and security.</p>
    <![endif]-->
    <header class="banner container">
        <h1 class="title">Diachromatic: Preparation and Quality Control of Hi-C and Capture Hi-C Data</h1>
        <div class="nav" id="navi">
            <ul class="nav-icons-home">
                <li id="home"><i class="fas fa-home"></i></li>
            </ul>
            <ul class="nav-icons-support">
                <li id="about"><i class="fas fa-info-circle"></i></li>
                <li id="settings"><i class="fas fa-cogs"></i></li>
                <li id="help"><i class="fas fa-question-circle"></i></li>
            </ul>
        </div>
    </header>
    <div class="container main">
        <div class="section sample">
            <table class="info-table">
                <tr>
                    <td class="table-label">
                        Date
                    </td>
                    <td class="table-data">
                        Sunday March 30, 2019 7:30pm
                    </td>
                </tr>
                <tr>
                    <td class="table-label">
                        Input (forward)
                    </td>
                    <td class="table-data">
                        ${input_fastq1!"n/a"}
                    </td>
                </tr>
                <tr>
                    <td class="table-label">
                        Input (reverse)
                    </td>
                    <td class="table-data">
                         ${input_fastq2!"n/a"}
                    </td>
                </tr>
            </table>
        </div>
        <div class="section extra">
            <div class="about">
                <h2>About</h2>
                <p>Diachromatic is a tool for ....</p>
            </div>
            <div class="settings">
                <h2>Settings</h2>
                <p>Diachromatic is a tool for ....</p>
            </div>
            <div class="help">
                <h2>Help</h2>
                <p>Diachromatic is a tool for ....</p>
            </div>
        </div>
        <!-- Section on truncation -->
        <div class="section">
             <a name="truncation"></a>
             <h2>Truncation statistics</h2>
             <p>Reads were processed to recognize filled in restriction sites of the enzyme
                ${restriction_enzyme!"n/a"}, i.e., ${filled_end_sequence!"n/a"}.</p>
             <p>Truncation was performed with a length threshold of ${length_threshold!"n/a"} nucleotides (reads whose length is below this threshold
                after truncation are removed, denoted <i>Too short to map</i> in the table). A total of ${removed_pairs_one_or_two_reads_too_short!"n/a"}
                 pairs contained either a forward or a reverse read (or both) that were too short to map.
             </p>
	     <br>
             <table class="redTable">
               <tr>
			<th></th><th>Forward Read</th><th>Reverse Read</th></tr>
               <tr>
                   <td><b>Total Reads</b></td>
                   <td>${total_read_pairs_processed!"n/a"}</td>
                   <td>${total_read_pairs_processed!"n/a"}</td>
               </tr>
               <tr>
                   <td><b>Truncated Reads</b></td>
                   <td>${truncated_forward_reads!"n/a"}</td>
                   <td>${truncated_reverse_reads!"n/a"}</td>
               </tr>
               <tr>
                   <td><b>Dangling Reads</b></td>
                   <td>${dangling_forward_reads!"n/a"}</td>
                   <td>${dangling_reverse_reads!"n/a"}</td>
               </tr>
               <tr>
                   <td><b>Too short to map</b></td>
                   <td>${short_removed_forward_reads!"n/a"}</td>
                   <td>${short_removed_reverse_reads!"n/a"}</td>
               </tr>
             </table>
        </div>
	<!-- Section on truncation - end -->
        <!-- Section on alignment -->
	<div class="section">
	<h2>Alignment statistics</h2>
		<p>
                Diachromatic performs two independent alignments of the truncated R1 and R2 reads
                with <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="__blank">bowtie2</a>.
                Read pairs for which either R1 or R2 cannot be mapped or cannot be mapped to a unique location are
                discarded. The remaining read pairs are re-paired and subjected to deduplication,
                i.e. for given read pairs with identical orientation and identically mapped 5' end positions
                only one read is used for downstream analysis.
		</p>
		<br>
		<table class="redTable">
			<tr>
			    <th></th><th>R1</th><th>R2</th><th>Pairs</th>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Total</b></td>
			    <td></td>
			    <td></td>
			    <td>${align_total_read_pairs_processed!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Unmapped</b></td>
			    <td>${align_unmapped_R1_reads!"n/a"}</td>
			    <td>${align_unmapped_R2_reads!"n/a"}</td>
			    <td>${align_unmapped_read_pairs!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Multimapped</b></td>
			    <td>${align_multimapped_R1_reads!"n/a"}</td>
			    <td>${align_multimapped_R2_reads!"n/a"}</td>
			    <td>${align_multimapped_read_pairs!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Paired</b></td>
			    <td></td>
			    <td></td>
			    <td>${align_paired_read_pairs!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Paired duplicated</b></td>
			    <td></td>
			    <td></td>
			    <td>${align_duplicated_pairs!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Paired unique</b></td>
			    <td></td>
			    <td></td>
			    <td>${align_unique_paired_read_pairs!"n/a"}</td>
			</tr>
		</table>
		<br>
                <div id="container_alignReadGraph" style="min-width: 150px; height: 350px; margin: 0 auto"></div>
	</div>
	<!-- Section on alignment - end -->
        <!-- Fragment statistics -->
        <div class="section">
	<h2>Fragment statistics</h2>
		<p>
                The paired unique pairs are subdivided into four disjoint categories: <i>un-ligated</i>, <i>self-ligated</i>,
                <i>chimeric</i> and <i>strange internal</i>. The <i>un-ligated</i> and <i>self-ligated</i> categories
                are further subdivided depending on whether the two reads were mapped to the same digest or
                whether the categorization was based on one of the corresponding size threshold.
                Reads that are neither <i>un-ligated</i> or <i>self-ligated</i> ligation are categorized as <i>chimeric</i>.
                The <i>chimeric</i> category id further subdivided into <i>too short</i>, <i>valid</i>, and <i>too long</i>
                depending on whether the calculated size is within the user defined range.
                The last category <i>strange internal</i> is a container for read pairs that cannot be explained by
                the processing logic. These pairs map to the same strand and, therefore, cannot be <i>un-ligated</i>
                or <i>self-ligated</i> but they are not chimeric because both reads map to the same digest.
                For all datasets that were analyzed, <i>strange internal</i> constitute only a small fraction of
                all paired unique pairs.
		</p>
		<br>
		<table class="redTable">
			<tr>
			    <th>Fragment type</th><th>Count</th>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Un-ligated</b></td>
			    <td>${align_unligated!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Un-ligated by size</b></td>
			    <td>${align_unligated_by_size!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"> Un-ligated by same internal</td>
			    <td>${align_unligated_same_internal!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Self-ligated</b></td>
			    <td>${align_self_ligated!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Self-ligated by size</td>
			    <td>${align_self_ligated_by_size!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Self-ligated by same internal</td>
			    <td>${align_self_ligated_same_internal!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Chimeric</b></td>
			    <td>${align_chimeric!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Chimeric too short</td>
			    <td>${align_chimeric_short!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Chimeric too long</td>
			    <td>${align_chimeric_long!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left">Chimeric valid</td>
			    <td>${align_chimeric_valid!"n/a"}</td>
			</tr>
			<tr>
			    <td style="text-align:left"><b>Strange internal</b></td>
			    <td>${align_strange_internal!"n/a"}</td>
			</tr>
		</table>
		<br>
                <div id="container_fragmentTypeCounts" style="min-width: 150px; height: 350px; margin: 0 auto"></div>
                TODO: Verify claim about strange internals. Revise text.
            </div>
	<!-- Fragment statistics -end -->
        <!-- Quality metrics -->
	<div class="section">
	<h2>Quality metrics</h2>
		<p>
                Explain quality metrics.
		</p>
		<table class="redTable">
			<tr>
			    <th>Quality metrics</th><th></th>
			</tr
			<tr><td style="text-align:left"><b> Yield of valid pairs (YVP)</b></td><td>${align_YVP!"n/a"}</td></tr>
			<tr><td style="text-align:left"><b> Cross-ligation coefficient (CLC)  </b></td><td>${align_CLC!"n/a"}</td></tr>
			<tr><td style="text-align:left"><b>Re-ligation coefficient (RLC)   </b></td><td>${align_RLC!"n/a"}</td></tr>
			<tr><td style="text-align:left"><b>Hi-C pair duplication rate (HPDR)</b></td><td>${align_HPDR!"n/a"}</td></tr>
		</table>
	</div>
	<!-- Quality metrics - end -->
	<!-- Detailed results -->
	<div class="section">
	<h2>Dangling end and trans pair subcategories</h2>
		<p>
		Dangling ends and trans read pairs do not constitute a separate category but may occur in different categories. The following table and bar chart show the proportion of
		dangling end and trans read pairs within the individual categories (the category <i>Chimeric</i> is broken down into <i>Chimeric - Too short</i>, <i>Chimeric - Valid size</i>,
		and <i>Chimeric - Too long</i>). Note: Trans pairs cannot occur within the categories <i>Un-ligated</i>, <i>Self-ligated</i> and <i>Strange internal</i> by defiition.
            	</p>
		<br>
		<table class="redTable">
			<tr>
				<th>Fragment type</th><th>Dangling</th><th>Trans</th>
			</tr>
			<tr>
				<td style="text-align:left"><b>Un-ligated</b></td><td>${align_n_paired_unique_un_ligated_dangling!"n/a"}</td><td>${align_n_paired_unique_un_ligated_trans!"n/a"}</td>
			</tr>
			<tr>
				<td style="text-align:left"><b>Self-ligated</b></td><td>${align_n_paired_unique_self_ligated_dangling!"n/a"}<td>${align_n_paired_unique_self_ligated_trans!"n/a"}</td></td>
			</tr>
			<tr>
				<td style="text-align:left"><b>Chimeric - Too short</b></td><td>${align_n_paired_unique_too_short_dangling!"n/a"}</td><td>${align_n_paired_unique_too_short_trans!"n/a"}</td>
			</tr>
			<tr>
				<td style="text-align:left"><b>Chimeric - Valid size</b></td><td>${align_n_paired_unique_valid_dangling!"n/a"}</td><td>${align_n_paired_unique_valid_trans!"n/a"}</td>
			</tr>
			<tr>
				<td style="text-align:left"><b>Chimeric - Too long</b></td><td>${align_n_paired_unique_too_long_dangling!"n/a"}</td><td>${align_n_paired_unique_too_long_trans!"n/a"}</td>
			</tr>
			<tr>
				<td style="text-align:left"><b>Strange internal</b></td><td>${align_n_paired_strange_internal_dangling!"n/a"}</td><td>${align_n_paired_strange_internal_trans!"n/a"}</td>
			</tr>
		</table>
		<br>
		<div id="container_danglingTransSubcategoryCounts" style="min-width: 150px; height: 350px; margin: 0 auto"></div>
	</div>
	<!-- Detailed results - end -->
	<!-- Fragment sizes -->
	<div class="section">
		<div id="container_fragmentSizeCounts" style="min-width: 150px; height: 350px; margin: 0 auto"></div>
	</div>
        <#if count??>
            <section>
                <article>
                    <a name="count"></a>
                    <h2>Count statistics</h2>
                    <ol>
                        <#list count as line>
                            <li>${line}</li>
                        </#list>
                    </ol>


                </article>
            </section>
        </#if> <!-- count section -->
    </div>
    <footer class="container">
      <p>Diachromatic &copy; 2019</p>
    </footer>
    <script type="text/javascript">
        <#include "jquery.js">
    </script>
    <script>
        <#include "highcharts.js">
    </script>
    <script type="text/javascript">
        <#include "custom-charts.js">
    </script>
    </body>
</html>

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
                        Sample
                    </td>
                    <td class="table-data">
                        Z17X
                    </td>
                </tr>
                <tr>
                    <td class="table-label">
                        Input Files
                    </td>
                    <td class="table-data">
                        z17x.align.stats.txt, z17x.truncation.stats.txt
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
             <h2>Truncation Statistics</h2>
             <p>Truncation was performed with a length threshold of ${length_threshold!"n/a"} nucleotides (reads whose length is below this threshold
                after truncation are removed, denoted <i>Too short to map</i> in the table). A total of ${removed_pairs_one_or_two_reads_too_short!"n/a"}
                 pairs contained either a forward or a reverse read (or both) that were too short to map.
             </p>
             <table class="redTable">
               <tr><th></th><th>Forward Read</th><th>Reverse Read</th></tr>
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
            <div id="container" style="min-width: 310px; height: 400px; margin: 0 auto"></div>
            <span id="align"></span>
        </div>
         <!-- Section on alignment -->
        <div class="section">
            <h2>Alignment Statistics</h2>
            <p>
                Diachromatic performs alignment of the truncated reads
                with <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="__blank">bowtie2</a>. TODO-more text.
            </p>
            <!-- Read Pair Section -->
            <div class="row">
                <div class="col-lg-5">
                    <h3>Read pair analysis</h3>
                    <p>TODO Some text..</p>
                    <table class="redTable">
                        <tr>
                            <th></th><th>Count</th>
                        </tr>
                        <tr>
                            <td><b>Total Reads</b></td>
                            <td>${align_total_read_pairs_processed!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Unmapped Read pairs</b></td>
                            <td>${align_unmapped_read_pairs!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Multimapped Read Pairs</b></td>
                            <td>${align_multimapped_read_pairs!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Paired Read pairs</b></td>
                            <td>${align_paired_read_pairs!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Unique paired Read pairs</b></td>
                            <td>${align_unique_paired_read_pairs!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Duplicated Read pairss</b></td>
                            <td>${align_duplicated_pairs!"n/a"}</td>
                        </tr>
                    </table>
                </div>
                <div class="col-lg-6">
                    <div id="container_alignReadPairGraph" style="min-width: 210px; height: 200px; margin: 0 auto"></div>
                </div>
            </div>
            <hr/>
            <div class="row">
                <!-- Read Section -->
                <div class="col-lg-5">
                    <h3>Read analysis</h3>
                    <p>TODO Some text about the reads..</p>
                    <table class="redTable">
                        <tr>
                            <th></th><th>Forward Read</th><th>Reverse Read</th>
                        </tr>
                        <tr>
                            <td><b>Unmapped Reads</b></td>
                            <td>${align_unmapped_R1_reads!"n/a"}</td>
                            <td>${align_unmapped_R2_reads!"n/a"}</td>
                        </tr>
                        <tr>
                            <td><b>Multimapped Reads</b></td>
                            <td>${align_multimapped_R1_reads!"n/a"}</td>
                            <td>${align_multimapped_R2_reads!"n/a"}</td>
                        </tr>
                    </table>
                </div>
                <div class="col-lg-6">
                    <div id="container_alignReadGraph" style="min-width: 150px; height: 200px; margin: 0 auto"></div>
                </div>
            </div>
            <span id="artefact"></span>
           <p> TODO -- consider reformating the table because many of the items are for read-pairs and not reads</p>
        </div>
        <!-- Artefact statistics -->
        <div class="section">
            <h2>Artefact Statistics</h2>
            <p>
                Artefacts are TODO--Write a summary here
            </p>
            <table class="redTable">
                <tr>
                    <th>Artefact type</th><th>Count</th>
                </tr>
                <tr>
                    <td><b>Unligated</b></td>
                    <td>${align_unligated!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>unligated_by_size</b></td>
                    <td>${align_unligated_by_size!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>unligated (same_internal)</b></td>
                    <td>${align_unligated_same_internal!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>self-ligated</b></td>
                    <td>${align_self_ligated!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>self-ligated_by_size</b></td>
                    <td>${align_self_ligated_by_size!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>self-ligated (same internal)</b></td>
                    <td>${align_self_ligated_same_internal!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>chimeric</b></td>
                    <td>${align_chimeric!"n/a"}</td>
                </tr>
                <tr>
                    <td><b> chimeric short</b></td>
                    <td>${align_chimeric_short!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>chimeric long</b></td>
                    <td>${align_chimeric_long!"n/a"}</td>
                </tr>
                <tr>
                    <td><b>chimeric valid</b></td>
                    <td>${align_chimeric_valid!"n/a"}</td>
                </tr>
                <tr>
                    <td><b> strange internal </b></td>
                    <td>${align_strange_internal!"n/a"}</td>
                </tr>
                <tr>
                    <td><b> dangling_end_pairs_total  </b></td>
                    <td>${align_dangling_end_pairs_total!"n/a"}</td>
                </tr>
                <tr>
                    <td><b> trans_pairs_total  </b></td>
                    <td>${align_trans_pairs_total!"n/a"}</td>
                </tr>
            </table>
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
            <br/><br/>
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
            <br/><br/>
            <table class="redTable">
                <caption>Summary</caption>
                <tr><th>Artefact type</th><th>Count</th></tr>
                <tr><td><b> Yield of valid pairs (YVP)</b></td><td>${align_YVP!"n/a"}</td></tr>
                <tr><td><b> Cross-ligation coefficient (CLC)  </b></td><td>${align_CLC!"n/a"}</td></tr>
                <tr><td><b>Re-ligation coefficient (RLC)   </b></td><td>${align_RLC!"n/a"}</td></tr>
                <tr><td><b>Hi-C pair duplication rate (HPDR)</b></td><td>${align_HPDR!"n/a"}</td></tr>
            </table>
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

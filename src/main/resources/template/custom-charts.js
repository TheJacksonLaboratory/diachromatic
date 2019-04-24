// Replace all
String.prototype.replaceAll = function(search, replacement) {
    var target = this;
    return target.split(search).join(replacement);
};
/* Method to convert array of strings to array of integers */
var convertToIntArray = function(ray){
    return ray.map(function (x) {
        x = x.replaceAll(",", "");
        return parseInt(x, 10);
    });
};

/* Alignment Charts */
var forwardReadCounts = ['${align_total_read_pairs_processed}','${unmapped_R1_reads}', '${multimapped_R1_reads}', '${align_paired_read_pairs}', '${align_duplicated_pairs}', '${align_unique_paired_read_pairs}'];
var reverseReadCounts = ['${align_total_read_pairs_processed}','${unmapped_R2_reads}', '${multimapped_R2_reads}', '${align_paired_read_pairs}', '${align_duplicated_pairs}', '${align_unique_paired_read_pairs}'];
var pairCounts = ['${align_total_read_pairs_processed}','${unmapped_read_pairs}', '${multimapped_read_pairs}', '${align_paired_read_pairs}', '${align_duplicated_pairs}', '${align_unique_paired_read_pairs}'];
Highcharts.chart('container_alignReadGraph', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Alignment'
    },
    subtitle: {
        text: 'Read and read pair counts following alignment'
    },
    xAxis: {
        categories: ['Total', 'Unmapped', 'Multimapped', 'Paired', 'Paired duplicated', 'Paired unique'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of reads/pairs'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [
        {
            name: 'R1',
            color: '#e69696',
            data: convertToIntArray(forwardReadCounts)
        },
        {
            name: 'R2',
            color: '#9696e6',
            data: convertToIntArray(reverseReadCounts)
        },
        {
            name: 'Pairs',
            color: '#743562',
            data: convertToIntArray(pairCounts)
        }
    ]
});

/* Interaction count statistics */
var interactionCounts = [
    '${count_total_interacting_fragments}',
    '${count_selected_interacting_fragments}',
    '${count_total_interaction_count}',
    '${count_interactions_between_selected_fragments}',
    '${count_interactions_between_selected_and_unselected_fragments}',
    '${count_interactions_between_unselected_fragments}'
    ];
Highcharts.chart('container_interactionCountsBarChart', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Interaction counts'
    },
    subtitle: {
        text: 'The total number of interactions (black) and numbers for the different subcategories (orange/gray)'
    },
    xAxis: {
        categories: ['Interacting digests', 'Active', 'Interactions', 'Active/Active',  'Active/Inactive', 'Inactive/Inactive'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of interactions'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        series: {
            colorByPoint: true
        },
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [
    {
        showInLegend: false,
        name: 'Digest and interaction counts',
        data: convertToIntArray(interactionCounts),
        colors: ['black', '#ffa500', 'black', '#ffa500', '#d1a452', '#a3a3a3']
    }
    ]
});




var transInteractionCounts = ['${count_n_singleton_interactions_trans}', '${count_n_gt1_interaction_count_trans}'];
var shortInteractionCounts = ['${count_n_singleton_interactions_short_range}', '${count_n_gt1_interaction_count_short_range}'];
var longInteractionCounts = ['${count_n_singleton_interactions_long_range}', '${count_n_gt1_interaction_count_long_range}'];
Highcharts.chart('container_singletonInteractions', {
  chart: {
    type: 'bar'
  },
  title: {
    text: 'Singleton interactions vs. others'
  },
  xAxis: {
    categories: ['Singleton', 'Other']
  },
  yAxis: {
    min: 0,
    title: {
      text: 'Composition of singleton interactions and other interactions'
    }
  },
  legend: {
    reversed: true
  },
  plotOptions: {
    series: {
      stacking: 'normal'
    }
  },
  series: [{
    name: 'Trans',
    data: convertToIntArray(transInteractionCounts)
  }, {
    name: 'Cis, Short range',
    data: convertToIntArray(shortInteractionCounts)
  }, {
    name: 'Cis, Long range',
    data: convertToIntArray(longInteractionCounts)
  }]
});


/* Dangling subcategories */
var danglingSubcategoryCounts = [
    '${align_dangling_end_pairs_total}',
    '${align_n_paired_unique_un_ligated_dangling}',
    '${align_n_paired_unique_self_ligated_dangling}',
    '${align_n_paired_unique_too_short_dangling}',
    '${align_n_paired_unique_valid_dangling}',
    '${align_n_paired_unique_too_long_dangling}',
    '${align_n_paired_strange_internal_dangling}'
    ];
var transSubcategoryCounts = [
    '${align_trans_pairs_total}',
    '${align_n_paired_unique_un_ligated_trans}',
    '${align_n_paired_unique_self_ligated_trans}',
    '${align_n_paired_unique_too_short_trans}',
    '${align_n_paired_unique_valid_trans}',
    '${align_n_paired_unique_too_long_trans}',
    '${align_n_paired_strange_internal_trans}'
    ];
var detailedCategoryCounts = [
    '${align_unique_paired_read_pairs}',
    '${align_unligated}',
    '${align_self_ligated}',
    '${align_chimeric_short}',
    '${align_chimeric_valid}',
    '${align_chimeric_long}',
    '${align_strange_internal}'
    ];
Highcharts.chart('container_danglingTransSubcategoryCounts', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Dangling and trans subcategories'
    },
    subtitle: {
        text: 'Fractions of dangling end and trans pairs'
    },
    xAxis: {
        categories: ['Total', 'Un-ligated', 'Self-ligated', 'Chimeric - Too short', 'Chimeric - Valid size', 'Chimeric - Too long', 'Strange internal'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of pairs'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [
        {
            name: 'Total',
            color: '#743562',
            data: convertToIntArray(detailedCategoryCounts)
        },
        {
            name: 'Dangling',
            color: '#00cc99',
            data: convertToIntArray(danglingSubcategoryCounts)
        },
        {
            name: 'Trans',
            color: '#3366cc',
            data: convertToIntArray(transSubcategoryCounts)
        }
    ]
});

/* Fragments size plots */

var ChimericFragmentSizeCounts = [${align_chimeric_fragment_size_count_array}];
var ChimericFragmentSizeCountsActive = [${align_chimeric_fragment_size_active_count_array}];
var UnLigatedFragmentSizeCountsActive = [${align_un_ligated_fragment_size_count_array}];
Highcharts.chart('container_fragmentSizeCounts', {
    title: {
        text: 'Distribution of fragment sizes'
    },
    chart: {
        type: 'line'
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Fragment counts'
        }
    },
    series: [
    {
        name: 'Chimeric',
        data: ChimericFragmentSizeCounts,
        pointStart: 1,
        color: '#9ADAFF',
        lineWidth: 0.5
    },
    {
        name: 'Chimeric active',
        data: ChimericFragmentSizeCountsActive,
        color: '#FFA500',
        pointStart: 1,
        lineWidth: 0.5
    },
    {
        name: 'Un-ligated',
        data: UnLigatedFragmentSizeCountsActive,
        color: '#A3A3A3',
        pointStart: 1,
        lineWidth: 0.5
    }
    ]
});

var SelfLigatedDigestSizeCounts = [${align_self_ligated_fragment_size_count_array}];
Highcharts.chart('container_SelfLigatedDigestSizeCounts', {
    title: {
        text: 'Distribution of self-ligating digest sizes'
    },
    chart: {
        type: 'line'
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Digest counts'
        }
    },
    series: [
    {
        name: 'Self-ligated',
        data: SelfLigatedDigestSizeCounts,
        pointStart: 1,
        color: '#9ADAFF',
        lineWidth: 0.5
    }
    ]
});

/* Fragment type count chart */
var fragmentTypeCounts = ['${align_unique_paired_read_pairs}', '${align_unligated}','${align_self_ligated}', '${align_chimeric}', '${align_strange_internal}'];
Highcharts.chart('container_fragmentTypeCounts', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Fragment types'
    },
    subtitle: {
        text: 'Read pair counts for different fragment types'
    },
    xAxis: {
        categories: ['Total', 'Un-ligated',  'Self-ligated', 'Chimeric', 'Strange internal'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of pairs'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [
    {
        name: 'Pairs',
        data: convertToIntArray(fragmentTypeCounts),
        color: '#743562'
    }
    ]
});

/* Trans/Cis scatterplot */
var transCisScatterValues = ${align_trans_cis_scatter_values_array};
var globalCLC = parseFloat(${align_global_clc});
Highcharts.chart('container_transCisScatterPlot', {
    chart: {
        type: 'scatter',
        zoomType: 'xy'
    },
    title: {
        text: 'Number of digest vs. proportion of trans reads for each chromosome'
    },
    plotOptions: {
        scatter: {
            dataLabels: {
                format: "{point.name}",
                allowOverlap: true,
                enabled: true
            }
        }
    },
    xAxis: {
        title: {
            text: 'Proportion of trans reads'
        },
        plotLines: [{
            color: 'gray', // Color value
            dashStyle: 'Dash',
            value: ${align_global_clc},
            width: 1,
        }]
    },
    yAxis: {
        title: {
            text: 'Number of digest per chromosome'
        }
    },
    series: [{
        regression: true,
        regressionSettings: {
        type: 'linear',
        color: 'red',
        dashStyle: 'Solid'
        },
        name: "Chromosome",
        color: '#9ADAFF',
        data: transCisScatterValues
    }

]
});



var readPairData = ['${total_read_pairs_processed}', '${unmapped_read_pairs}', '${multimapped_read_pairs}',
                    '${paired_read_pairs}', '${unique_paired_read_pairs}', '${duplicated_pairs}'];
Highcharts.chart('container_alignReadPairGraph', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Alignment'
    },
    subtitle: {
        text: 'Read and pair counts following alignment'
    },
    xAxis: {
        categories: [  'Total readpairs',  'Unmapped readpairs', 'Multimapped readpairs',
            'Paired readpairs', 'Unique paired read pairs', 'Duplicated pairs'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of reads/pairs'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [{
        name: 'Count',
        data: convertToIntArray(readPairData)
    }]
});

var truncationForwardData = ['${total_read_pairs_processed}', '${truncated_forward_reads}', '${dangling_forward_reads}',
                        '${short_removed_forward_reads}'];
var truncationReverseData = ['${total_read_pairs_processed}', '${truncated_reverse_reads}', '${dangling_reverse_reads}',
                        '${short_removed_reverse_reads}'];




/* Truncation Charts */
Highcharts.chart('container', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Truncation'
    },
    subtitle: {
        text: 'read counts following truncation'
    },
    xAxis: {
        categories: ['Total reads', 'Truncated', 'Dangling', 'Too short'],
        crosshair: true
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Number of reads'
        }
    },
    tooltip: {
        headerFormat: '<span style="font-size:10px">{point.key}</span><table>',
        pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
            '<td style="padding:0"><b>{point.y:.1f} mm</b></td></tr>',
        footerFormat: '</table>',
        shared: true,
        useHTML: true
    },
    plotOptions: {
        column: {
            pointPadding: 0.2,
            borderWidth: 0
        }
    },
    series: [{
        name: 'Forward read',
        data: convertToIntArray(truncationForwardData)
    }, {
        name: 'Reverse read',
        data: convertToIntArray(truncationReverseData)
    }]
});

$(".nav-icons-support li").click(function(event) {
    var target = event.currentTarget.id;
    reset();
    $(".section").not(".extra").hide();
    $(".section.extra").show();
    console.log(event);
    if(target === "about"){
        $(".about").show();
    } else if(target ===  "settings"){
        $(".settings").show();
    } else if(target === "help"){
        $(".help").show();
    }

});

$(".nav-icons-home li").click(function(){
    $(".section.extra").hide();
    $(".section").not(".extra").show();
});

var reset = function(){
    $(".about").hide();
    $(".settings").hide();
    $(".help").hide();
};


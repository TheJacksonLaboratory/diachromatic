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

/* Artifact Charts */
var artifactDataSimple = ['${align_unique_paired_read_pairs}', '${align_unligated}','${align_self_ligated}', '${align_chimeric}', '${align_strange_internal}'];
Highcharts.chart('container_artifactCountsSimpleGraph', {
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
        data: convertToIntArray(artifactDataSimple),
        color: '#743562'
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


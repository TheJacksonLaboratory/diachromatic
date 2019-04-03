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
var forwardReadCounts = ['${unmapped_R1_reads}', '${multimapped_R1_reads}'];
var reverseReadCounts = ['${unmapped_R2_reads}', '${multimapped_R2_reads}'];
Highcharts.chart('container_alignReadGraph', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Read Alignment'
    },
    subtitle: {
        text: 'read counts following alignment'
    },
    xAxis: {
        categories: ['Unmapped reads', 'Multimapped reads'],
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
    series: [
        {
            name: 'Forward read',
            data: convertToIntArray(forwardReadCounts)
        },
        {
            name: 'Reverse read',
            data: convertToIntArray(reverseReadCounts)
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
        text: 'Read Pair Alignment'
    },
    subtitle: {
        text: 'read counts following alignment'
    },
    xAxis: {
        categories: [  'Total readpairs',  'Unmapped readpairs', 'Multimapped readpairs',
            'Paired readpairs', 'Unique paired read pairs', 'Duplicated pairs'],
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
        name: 'Count',
        data: convertToIntArray(readPairData)
    }]
});

var truncationForwardData = ['${total_read_pairs_processed}', '${truncated_forward_reads}', '${dangling_forward_reads}',
                        '${short_removed_forward_reads}'];
var truncationReverseData = ['${total_read_pairs_processed}', '${truncated_reverse_reads}', '${dangling_reverse_reads}',
                        '${short_removed_reverse_reads}'];
/** Truncation Alignment */
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


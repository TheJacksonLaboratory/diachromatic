package org.jax.diachromatic.summarize;

class TruncationJavaScript {
    private final int total_read_pairs_processed;
    private final int truncated_forward_reads;
    private final int truncated_reverse_reads;
    private final int dangling_forward_reads;
    private final int dangling_reverse_reads;
    private final int short_removed_forward_reads;
    private final int short_removed_reverse_reads;

    TruncationJavaScript(int readpairs,int truncatedF, int truncatedR, int danglingF, int danglingR,
                         int shortF, int shortR) {
        this.total_read_pairs_processed=readpairs;
        this.truncated_forward_reads=truncatedF;
        this.truncated_reverse_reads=truncatedR;
        this.dangling_forward_reads=danglingF;
        this.dangling_reverse_reads=danglingR;
        this.short_removed_forward_reads=shortF;
        this.short_removed_reverse_reads=shortR;
    }


    public String getJavaScript() {
        StringBuilder sb = new StringBuilder();
        sb.append("Highcharts.chart('container', {\n" +
                "    chart: {\n" +
                "        type: 'column'\n" +
                "    },\n" +
                "    title: {\n" +
                "        text: 'Truncation'\n" +
                "    },\n" +
                "    subtitle: {\n" +
                "        text: 'read counts following truncation'\n" +
                "    },");
        sb.append(" xAxis: {\n" +
                "        categories: [ " +
                "            'Total reads', " +
                "            'Not truncated'," +
                "            'Truncated'," +
                "            'Too short'],\n" +
                "crosshair: true " +
                "},\n");
        sb.append(" yAxis: {\n" +
                "        min: 0,\n" +
                "        title: {\n" +
                "            text: 'Number of reads'\n" +
                "        }\n" +
                "    },");
        sb.append("tooltip: {\n" +
                "        headerFormat: '<span style=\"font-size:10px\">{point.key}</span><table>',\n" +
                "        pointFormat: '<tr><td style=\"color:{series.color};padding:0\">{series.name}: </td>' +\n" +
                "            '<td style=\"padding:0\"><b>{point.y:.1f} mm</b></td></tr>',\n" +
                "        footerFormat: '</table>',\n" +
                "        shared: true,\n" +
                "        useHTML: true\n" +
                "    },\n");
        sb.append(" plotOptions: {\n" +
                "        column: {\n" +
                "            pointPadding: 0.2,\n" +
                "            borderWidth: 0\n" +
                "        }\n" +
                "    },");

        sb.append(String.format("series: [{\n" +
                "        name: 'Forward read',\n" +
                "        data: [%d, %d,%d,%d]\n" +
                "    }, ",total_read_pairs_processed,truncated_forward_reads,dangling_forward_reads,short_removed_forward_reads));
        sb.append(String.format("{\n" +
                "        name: 'Reverse read',\n" +
                "        data: [%d, %d,%d,%d]\n" +
                "     }]\n" +
                "});",total_read_pairs_processed,truncated_reverse_reads,dangling_reverse_reads,short_removed_reverse_reads));

        return sb.toString();
    }

}

/*
     name: 'Berlin',
        data: [42.4, 33.2, 34.5, 39.7, 52.6, 75.5, 57.4, 60.4, 47.6, 39.1, 46.8, 51.1]

    }]
});



	// Mapping Analysis
		$(function () {
			$('#truncating_mapping_plot').highcharts({
				colors: [
				   '#0d233a',
				   '#2f7ed8',
				   '#8bbc21',
				   '#910000',
				   '#1aadce',
				   '#492970',
				   '#f28f43',
				   '#77a1e5',
				   '#c42525',
				   '#a6c96a'
				],
				chart: {
					type: 'column',
					marginRight:0
				},
				title: {
					text: ''
				},
				xAxis: {
					categories: ['Total Reads',
				                                 'Not Truncated',
					                         'Truncated',
                                                                 'Too short to map',
								 'Unique Alignments',
								 'Multiple Alignments',
								 'Failed To Align',
								 'Paired']
				},
				yAxis: {
					title: {
						text: 'Number of Reads'
					}
				},
				credits: {
					enabled: false
				},
				tooltip: {
					formatter: function() {
						if(this.series.name == 'Read 1'){
							return '<b>'+ this.series.name + '</b><br/>'+
								this.x + ': ' +((this.y / 99742) * 100).toFixed(1) + '%';
						} else {
							return '<b>'+ this.series.name + '</b><br/>'+
								this.x + ': ' +((this.y / 99742) * 100).toFixed(1) + '%';
						}
					}
				},
				legend: {
					align: 'right',
					verticalAlign: 'top',
					y: 10
				},
				series: [{
					name: 'Read 1',
					data: [99742, 98802, 940, 333, 71235, 16631, 11543, 54021]
				}, {
					name: 'Read 2',
					data: [99742, 98875, 867, 339, 69787, 17419, 12197, 54021]
				}]
			});
		});
 */

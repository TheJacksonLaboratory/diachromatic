package org.jax.diachromatic.summarize;

/**
 * Convenience class to create a high-chart compatible JavaScript function to display
 * the align results as a bar chart.
 */
class AlignJavaScript {
    final private int total_read_pairs_processed;
    final private int unmapped_read_pairs;
    final private int unmapped_R1_reads;
    final private int unmapped_R2_reads;
    final private int multimapped_read_pairs;
    final private int multimapped_R1_reads;
    final private int multimapped_R2_reads;
    final private int paired_read_pairs;
    final private int unique_paired_read_pairs;
    final private int duplicated_pairs;

    AlignJavaScript( int total_read_pairs_processed,
            int unmapped_read_pairs,
            int unmapped_R1_reads,
            int unmapped_R2_reads,
            int multimapped_read_pairs,
            int multimapped_R1_reads,
            int multimapped_R2_reads,
            int paired_read_pairs,
            int unique_paired_read_pairs,
            int duplicated_pairs) {
        this.total_read_pairs_processed=total_read_pairs_processed;
        this.unmapped_read_pairs=unmapped_read_pairs;
        this.unmapped_R1_reads=unmapped_R1_reads;
        this.unmapped_R2_reads=unmapped_R2_reads;
        this.multimapped_read_pairs=multimapped_read_pairs;
        this.multimapped_R1_reads=multimapped_R1_reads;
        this.multimapped_R2_reads=multimapped_R2_reads;
        this.paired_read_pairs=paired_read_pairs;
        this.unique_paired_read_pairs=unique_paired_read_pairs;
        this.duplicated_pairs=duplicated_pairs;
    }


    /**
     * Make a highchart script to show the statistics for read-pairs (without the data on reads)
     * @return string representing a JavaScript function to build a chart.
     */
    String getReadJavaScript() {
        StringBuilder sb = new StringBuilder();
        sb.append("Highcharts.chart('container_alignRead', {\n" +
                "    chart: {\n" +
                "        type: 'column'\n" +
                "    },\n" +
                "    title: {\n" +
                "        text: 'Alignment'\n" +
                "    },\n" +
                "    subtitle: {\n" +
                "        text: 'read counts following alignment'\n" +
                "    },width : { 500px },");
        sb.append(" xAxis: {\n" +
                "        categories: ['Unmapped reads', 'Multimapped reads'], "+
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
                        "        data: [%d, %d]\n" +
                        "    }, ",unmapped_R1_reads,multimapped_R1_reads));
        sb.append(String.format("{\n" +
                        "        name: 'Reverse read',\n" +
                        "        data: [%d, %d]\n" +
                        "     }]\n" +
                        "});",unmapped_R2_reads,multimapped_R2_reads));

        return sb.toString();
    }


    String getReadPairJavaScript() {
        StringBuilder sb = new StringBuilder();
        sb.append("Highcharts.chart('container_alignReadPair', {\n" +
                "    chart: {\n" +
                "        type: 'column'\n" +
                "    },\n" +
                "    title: {\n" +
                "        text: 'Alignment'\n" +
                "    },\n" +
                "    subtitle: {\n" +
                "        text: 'read counts following alignment'\n" +
                "    },");
        sb.append(" width: 400,");
        sb.append(" xAxis: {\n" +
                "        categories: [ " +
                " 'Total readpairs', " +
                " 'Unmapped readpairs'," +
                " 'Multimapped readpairs', \n" +
                " 'Paired readpairs'," +
                " 'Unique paired read pairs',"+
                " 'Duplicated pairs'], "+
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
                        "        name: 'Count',\n" +
                        "        data: [%d, %d,%d,%d,%d, %d]\n" +
                        "    }); ",total_read_pairs_processed,unmapped_read_pairs,multimapped_read_pairs,
                paired_read_pairs,unique_paired_read_pairs,duplicated_pairs));

        return sb.toString();
    }


    String getJavaScript() {
        StringBuilder sb = new StringBuilder();
        sb.append("Highcharts.chart('container_align', {\n" +
                "    chart: {\n" +
                "        type: 'column'\n" +
                "    },\n" +
                "    title: {\n" +
                "        text: 'Alignment'\n" +
                "    },\n" +
                "    subtitle: {\n" +
                "        text: 'read counts following alignment'\n" +
                "    },");
        sb.append(" xAxis: {\n" +
                "        categories: [ " +
                " 'Total reads', " +
                " 'Unmapped readpairs'," +
                " 'Unmapped reads'," +
               " 'Multimapped readpairs', \n" +
                " 'Multimapped reads'," +
                " 'Paired readpairs'," +
                " 'Unique paired read pairs',"+
                " 'Duplicated pairs'], "+
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
                "        data: [%d, %d,%d,%d,%d, %d,%d,%d]\n" +
                "    }, ",total_read_pairs_processed,unmapped_read_pairs,unmapped_R1_reads,multimapped_read_pairs,multimapped_R1_reads,
                paired_read_pairs,unique_paired_read_pairs,duplicated_pairs));
        sb.append(String.format("{\n" +
                "        name: 'Reverse read',\n" +
                "        data: [%d, %d,%d,%d,%d, %d,%d,%d]\n" +
                "     }]\n" +
                "});",total_read_pairs_processed,unmapped_read_pairs,unmapped_R2_reads,multimapped_read_pairs,multimapped_R2_reads,
                paired_read_pairs,unique_paired_read_pairs,duplicated_pairs));



        return sb.toString();
    }
}

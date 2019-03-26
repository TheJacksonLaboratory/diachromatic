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





    String getJavaScript() {
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


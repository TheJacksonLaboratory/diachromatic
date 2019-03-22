package org.jax.diachromatic.summarize;

import freemarker.template.Configuration;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import freemarker.template.Version;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Summarizer {
    private static final Logger logger = LogManager.getLogger();

    private final String truncatePath;
    private final String alignPath;
    private final String countPath;

    /** Map of data that will be used for the FreeMark template. */
    private final Map<String, Object> templateData= new HashMap<>();
    /** FreeMarker configuration object. */
    private final Configuration cfg;


    public Summarizer(String truncFile, String alignFile, String countFile) {
        this.truncatePath=truncFile;
        this.alignPath=alignFile;
        this.countPath=countFile;
        this.cfg = new Configuration(new Version("2.3.23"));
        cfg.setDefaultEncoding("UTF-8");
        ClassLoader classLoader = Summarizer.class.getClassLoader();
        cfg.setClassLoaderForTemplateLoading(classLoader,"");
        try {
            ingestJavaScriptFiles();
        } catch (IOException e) {
            logger.error("Could not ingest java script files: {}", e.getMessage());
            // nothing we can do. The HTML pages will not display the charts
            // this should happen very rarely to never though
        }
        if(truncFile!=null) {
            parseTruncateData();
        }
        if(alignFile!=null) {
            parseAlignData();
        }
        if(countFile!=null) {
            parseCountData();
        }
    }

    /**
     * Three javascript files for highcharts are in the src/main/resources/template directory.
     * We need to copy these into the HTML page for the charts to work. We basically read them
     * into one string each and "paste" them into the template using the templateData map.
     */
    private void ingestJavaScriptFiles() throws IOException{
        ClassLoader classLoader = Summarizer.class.getClassLoader();
        URL url = classLoader.getResource("template/highcharts.js");
        if (url==null) {
            logger.error("template/highcharts.js not found");
            // this will mean that the charts cannot be displayed on the HTML page
            // however, this error should "never" happen. Just return if for some reason it does
            return;
        }
        String highchartspath = url.getFile();
        String highcharts =  new String ( Files.readAllBytes( Paths.get(highchartspath) ) );
        this.templateData.put("highcharts",highcharts);
        url = classLoader.getResource("template/exporting.js");
        if (url==null) {
            logger.error("template/exporting.js not found");
            // this will mean that the charts cannot be displayed on the HTML page
            // however, this error should "never" happen. Just return if for some reason it does
            return;
        }
        String exportingpath = url.getFile();
        String exporting =  new String ( Files.readAllBytes( Paths.get(exportingpath) ) );
        this.templateData.put("exporting",exporting);
        url = classLoader.getResource("template/export-data.js");
        if (url==null) {
            logger.error("template/export-data.js not found");
            // this will mean that the charts cannot be displayed on the HTML page
            // however, this error should "never" happen. Just return if for some reason it does
            return;
        }
        String exportdatapath = url.getFile();
        String exportdata =  new String ( Files.readAllBytes( Paths.get(exportdatapath) ) );
        this.templateData.put("exportdata",exportdata);

    }


    private int getIntegerValue(String f) {
        f=f.trim();
        int UNINITIALIZED=-1;
        Pattern pattern1 = Pattern.compile("\\d+");
        Pattern pattern2 = Pattern.compile("(\\d+)\\s+\\(\\d+\\.\\d+%\\)");
        Matcher matcher1=pattern1.matcher(f);
        if (matcher1.matches()) {
            return Integer.parseInt(f);
        }
        Matcher matcher2 = pattern2.matcher(f);
        if (matcher2.matches()) {
            String g1 = matcher2.group(1);
            return Integer.parseInt(g1);
        }
        // if we get here, we could not match -- probably we need to change the code that
        // outputs the data to the text files.
        return UNINITIALIZED;
    }

    private void parseTruncateData() {
        logger.trace("Parsing the truncation data at {}", truncatePath);
        int UNINITIALIZED=-1;
        int total_read_pairs_processed=UNINITIALIZED;
        int truncated_forward_reads=UNINITIALIZED;
        int truncated_reverse_reads=UNINITIALIZED;
        int dangling_forward_reads=UNINITIALIZED;
        int dangling_reverse_reads=UNINITIALIZED;
        int short_removed_forward_reads=UNINITIALIZED;
        int short_removed_reverse_reads=UNINITIALIZED;

        try (
            BufferedReader br = new BufferedReader(new FileReader(truncatePath))) {
            String line;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                String []fields=line.split(":");
                if (fields.length!=2) continue; // skip non key-value lines, they are comments
                templateData.put(fields[0],fields[1]);
                switch (fields[0]) {
                    case "total_read_pairs_processed" :
                        total_read_pairs_processed = getIntegerValue(fields[1]);
                        break;
                    case "truncated_forward_reads":
                        truncated_forward_reads=getIntegerValue(fields[1]);
                        break;
                    case "truncated_reverse_reads":
                        truncated_reverse_reads=getIntegerValue(fields[1]);
                        break;
                    case "dangling_forward_reads":
                        dangling_forward_reads=getIntegerValue(fields[1]);
                        break;
                    case "dangling_reverse_reads":
                        dangling_reverse_reads=getIntegerValue(fields[1]);
                        break;
                    case "short_removed_forward_reads":
                        short_removed_forward_reads=getIntegerValue(fields[1]);
                        break;
                    case "short_removed_reverse_reads":
                        short_removed_reverse_reads=getIntegerValue(fields[1]);
                        break;
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        TruncationJavaScript js = new  TruncationJavaScript(total_read_pairs_processed,
                truncated_forward_reads,truncated_reverse_reads,
                dangling_forward_reads,dangling_reverse_reads,
                short_removed_forward_reads,short_removed_reverse_reads);
        templateData.put("truncationjs",js.getJavaScript());
        System.err.println(js.getJavaScript());
    }

    private void parseAlignData() {
        logger.trace("Parsing the align data at {}", alignPath);
        logger.trace("Parsing the truncation data at {}", truncatePath);
        int UNINITIALIZED=-1;
        int total_read_pairs_processed=UNINITIALIZED;
        int unmapped_read_pairs=UNINITIALIZED;
        int unmapped_R1_reads=UNINITIALIZED;
        int unmapped_R2_reads=UNINITIALIZED;
        int multimapped_read_pairs=UNINITIALIZED;
        int multimapped_R1_reads=UNINITIALIZED;
        int multimapped_R2_reads=UNINITIALIZED;
        int paired_read_pairs=UNINITIALIZED;
        int unique_paired_read_pairs=UNINITIALIZED;
        int duplicated_pairs=UNINITIALIZED;

       String line;
        try (BufferedReader br = new BufferedReader(new FileReader(alignPath))) {
        while ((line=br.readLine())!=null) {
            System.out.println(line);
            String[] fields = line.split(":");
            if (fields.length != 2) continue; // skip non key-value lines, they are comments
            templateData.put(String.format("align_%s",fields[0].trim()), fields[1].trim());
            System.err.println( String.format("align_%s",fields[0]));
            switch (fields[0]) {
                case "total_read_pairs_processed":
                    total_read_pairs_processed = getIntegerValue(fields[1]);
                    break;
                case "unmapped_read_pairs":
                    unmapped_read_pairs=getIntegerValue(fields[1]);
                    break;
                case "unmapped_R1_reads":
                    unmapped_R1_reads=getIntegerValue(fields[1]);
                   break;
                case "unmapped_R2_reads":
                    unmapped_R2_reads=getIntegerValue(fields[1]);
                    break;
                case "multimapped_read_pairs":
                    multimapped_read_pairs=getIntegerValue(fields[1]);
                    break;
                case "multimapped_R1_reads":
                    multimapped_R1_reads=getIntegerValue(fields[1]);
                    break;
                case "multimapped_R2_reads":
                    multimapped_R2_reads=getIntegerValue(fields[1]);
                    break;
                case "paired_read_pairs":
                    paired_read_pairs=getIntegerValue(fields[1]);
                    break;
                case "unique_paired_read_pairs":
                    unique_paired_read_pairs=getIntegerValue(fields[1]);
                    break;
                case "duplicated_pairs":
                    duplicated_pairs=getIntegerValue(fields[1]);
                    break;
            }
        }

    } catch (IOException e) {
        e.printStackTrace();
    }

    /*
    int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
        int =UNINITIALIZED;
     */
        AlignJavaScript js = new  AlignJavaScript(total_read_pairs_processed,
                unmapped_read_pairs,unmapped_R1_reads,
                unmapped_R2_reads,multimapped_read_pairs,
                multimapped_R1_reads,multimapped_R2_reads,
                paired_read_pairs,unique_paired_read_pairs,duplicated_pairs);
        templateData.put("alignjs",js.getJavaScript());
        System.err.println(js.getJavaScript());
    }

    private void parseCountData() {
        logger.trace("Parsing the count data at {}", countPath);
        List<String> countMap = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(countPath))) {
            String line;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                countMap.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace("Putting a total of {} items into the align count", countMap.size());
        templateData.put("count",countMap);
    }


    public void outputFile(String prefix){
        String outname=String.format("%s.summary.stats.html",prefix);
        logger.trace("Writing HTML file to {}",outname);
        try (BufferedWriter out = new BufferedWriter(new FileWriter(outname))) {
            Template template = cfg.getTemplate("template/diachromatic-html-template.ftl");
            template.process(templateData, out);
        } catch (TemplateException | IOException te) {
            te.printStackTrace();
        }
    }
}

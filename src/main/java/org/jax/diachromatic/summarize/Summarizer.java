package org.jax.diachromatic.summarize;

import freemarker.template.Configuration;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import freemarker.template.Version;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
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

    // Map of data that will be used for the FreeMark template.
    private final Map<String, Object> templateData = new HashMap<>();

    // FreeMarker configuration object.
    private final Configuration cfg;

    public Summarizer(String truncFile, String alignFile, String countFile) {
        this.truncatePath = truncFile;
        this.alignPath = alignFile;
        this.countPath = countFile;
        this.cfg = new Configuration(new Version("2.3.23"));
        cfg.setDefaultEncoding("UTF-8");
        ClassLoader classLoader = Summarizer.class.getClassLoader();
        cfg.setClassLoaderForTemplateLoading(classLoader, "");

        if (truncFile != null) {
            parseTruncateData();
        }
        if (alignFile != null) {
            parseAlignData();
        }
        if (countFile != null) {
            parseCountData();
        }
    }

    private int getIntegerValue(String f) {
        f = f.trim();
        int UNINITIALIZED = -1;
        Pattern pattern1 = Pattern.compile("\\d+");
        Pattern pattern2 = Pattern.compile("(\\d+)\\s+\\(\\d+\\.\\d+%\\)");
        Matcher matcher1 = pattern1.matcher(f);
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
        int UNINITIALIZED = -1;
        String EMPTY_STRING = "";
        int total_read_pairs_processed = UNINITIALIZED;
        int truncated_forward_reads = UNINITIALIZED;
        int truncated_reverse_reads = UNINITIALIZED;
        int dangling_forward_reads = UNINITIALIZED;
        int dangling_reverse_reads = UNINITIALIZED;
        int short_removed_forward_reads = UNINITIALIZED;
        int short_removed_reverse_reads = UNINITIALIZED;
        String restriction_enzyme = EMPTY_STRING;
        String filled_end_sequence = EMPTY_STRING;

        logger.trace("Parsing the truncation data at {}", truncatePath);
        try (BufferedReader br = new BufferedReader(new FileReader(truncatePath))) {
            String line;
            while ((line=br.readLine())!=null) {
                logger.debug(line);
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
                    case "restriction_enzyme":
                        restriction_enzyme=fields[1];
                        break;
                    case "filled_end_sequence":
                        filled_end_sequence = fields[1];
                        break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Truncation Data
        templateData.put("total_read_pairs_processed", total_read_pairs_processed);
        templateData.put("truncated_forward_reads", truncated_forward_reads);
        templateData.put("truncated_reverse_reads", truncated_reverse_reads);
        templateData.put("dangling_forward_reads", dangling_forward_reads);
        templateData.put("dangling_reverse_reads", dangling_reverse_reads);
        templateData.put("short_removed_forward_reads", short_removed_forward_reads);
        templateData.put("short_removed_reverse_reads", short_removed_reverse_reads);
    }

    private void parseAlignData() {
        int UNINITIALIZED = -1;
        int total_read_pairs_processed = UNINITIALIZED;
        int unmapped_read_pairs = UNINITIALIZED;
        int unmapped_R1_reads = UNINITIALIZED;
        int unmapped_R2_reads = UNINITIALIZED;
        int multimapped_read_pairs = UNINITIALIZED;
        int multimapped_R1_reads = UNINITIALIZED;
        int multimapped_R2_reads = UNINITIALIZED;
        int paired_read_pairs = UNINITIALIZED;
        int unique_paired_read_pairs = UNINITIALIZED;
        int duplicated_pairs = UNINITIALIZED;
        String line;

        logger.trace("Parsing the align data at {}", alignPath);
        logger.trace("Parsing the truncation data at {}", truncatePath);
        try (BufferedReader br = new BufferedReader(new FileReader(alignPath))) {
            while ((line=br.readLine())!=null) {
                logger.debug(line);
                String[] fields = line.split(":");
                if (fields.length != 2) continue; // skip non key-value lines, they are comments
                String[] fields2 = fields[1].split(" ");
                if(!fields[0].contains("YVP") && !fields[0].contains("CLC") && !fields[0].contains("RLC") && !fields[0].contains("HPDR")) {
                    if(fields[0].contains("array")) {
                        templateData.put(String.format("align_%s", fields[0].trim()), fields[1].replace("%",":"));
                        logger.trace(fields[1].length());
                    } else {
                        templateData.put(String.format("align_%s", fields[0].trim()), getIntegerValue(fields2[0].trim()));
                    }
                } else {
                        templateData.put(String.format("align_%s", fields[0].trim()), fields2[0].trim());
                }
                if(fields[0].contains("global_clc")) {
                    templateData.put(String.format("align_%s", fields[0].trim()), fields2[0].trim());
                }



                logger.debug(String.format("align_%s %s",fields[0].trim(), fields2[0].trim()));
                logger.error(String.format("align_%s",fields[0]));
                switch (fields[0]) {
                    case "total_read_pairs_processed":
                        total_read_pairs_processed = getIntegerValue(fields[1]);
                        break;
                    case "unmapped_read_pairs":
                        unmapped_read_pairs = getIntegerValue(fields[1]);
                        break;
                    case "unmapped_R1_reads":
                        unmapped_R1_reads = getIntegerValue(fields[1]);
                        break;
                    case "unmapped_R2_reads":
                        unmapped_R2_reads = getIntegerValue(fields[1]);
                        break;
                    case "multimapped_read_pairs":
                        multimapped_read_pairs = getIntegerValue(fields[1]);
                        break;
                    case "multimapped_R1_reads":
                        multimapped_R1_reads = getIntegerValue(fields[1]);
                        break;
                    case "multimapped_R2_reads":
                        multimapped_R2_reads = getIntegerValue(fields[1]);
                        break;
                    case "paired_read_pairs":
                        paired_read_pairs = getIntegerValue(fields[1]);
                        break;
                    case "unique_paired_read_pairs":
                        unique_paired_read_pairs = getIntegerValue(fields[1]);
                        break;
                    case "duplicated_pairs":
                        duplicated_pairs = getIntegerValue(fields[1]);
                        break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Data for read alignment
        templateData.put("unmapped_R1_reads", unmapped_R1_reads);
        templateData.put("unmapped_R2_reads", unmapped_R2_reads);
        templateData.put("multimapped_R1_reads", multimapped_R1_reads);
        templateData.put("multimapped_R2_reads", multimapped_R2_reads);

        // Data for read pair alignment
        templateData.put("total_read_pairs_processed", total_read_pairs_processed);
        templateData.put("unmapped_read_pairs", unmapped_read_pairs);
        templateData.put("multimapped_read_pairs", multimapped_read_pairs);
        templateData.put("paired_read_pairs", paired_read_pairs);
        templateData.put("unique_paired_read_pairs", unique_paired_read_pairs);
        templateData.put("duplicated_pairs", duplicated_pairs);
    }

    private void parseCountData() {
        List<String> countMap = new ArrayList<>();

        logger.trace("Parsing the count data at {}", countPath);
        try (BufferedReader br = new BufferedReader(new FileReader(countPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                System.out.println(line);
                countMap.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace("Putting a total of {} items into the align count", countMap.size());
        templateData.put("count", countMap);
    }

    public void outputFile(String prefix){
        String outname = String.format("%s.summary.stats.html", prefix);
        logger.trace("Writing HTML file to {}", outname);
        try (BufferedWriter out = new BufferedWriter(new FileWriter(outname))) {
            Template template = cfg.getTemplate("template/diachromatic-html-template.ftl");
            template.process(templateData, out);
        } catch (TemplateException | IOException te) {
            te.printStackTrace();
        }
    }
}

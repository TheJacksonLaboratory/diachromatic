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

public class Summarizer {
    private static final Logger logger = LogManager.getLogger();

    private final String truncatePath;
    private final String alignPath;
    private final String countPath;

    /** Map of data that will be used for the FreeMark template. */
    protected final Map<String, Object> templateData= new HashMap<>();
    /** FreeMarker configuration object. */
    protected final Configuration cfg;


    public Summarizer(String truncFile, String alignFile, String countFile) {
        this.truncatePath=truncFile;
        this.alignPath=alignFile;
        this.countPath=countFile;
        this.cfg = new Configuration(new Version("2.3.23"));
        cfg.setDefaultEncoding("UTF-8");
        ClassLoader classLoader = Summarizer.class.getClassLoader();
        cfg.setClassLoaderForTemplateLoading(classLoader,"");
        if(truncFile!=null) {
            parseTruncateData();
        }
        if(alignFile!=null) {
            parseAlignData();
        }
        if(alignFile!=null) {
            parseCountData();
        }
    }

    private void parseTruncateData() {
        logger.trace("Parsing the truncation data at {}", truncatePath);
        List<String> truncateMap = new ArrayList<>();
        try (
            BufferedReader br = new BufferedReader(new FileReader(truncatePath))) {
            String line;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                truncateMap.add(line);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace("Putting a total of {} items into the truncate map", truncateMap.size());
        templateData.put("truncate",truncateMap);
    }

    private void parseAlignData() {
        logger.trace("Parsing the align data at {}", alignPath);
        List<String> alignMap = new ArrayList<>();
        try (
                BufferedReader br = new BufferedReader(new FileReader(alignPath))) {
            String line;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                alignMap.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace("Putting a total of {} items into the align map", alignMap.size());
        templateData.put("align",alignMap);
    }

    private void parseCountData() {
        logger.trace("Parsing the count data at {}", countPath);
        List<String> countMap = new ArrayList<>();
        try (
                BufferedReader br = new BufferedReader(new FileReader(countPath))) {
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

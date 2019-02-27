package org.jax.diachromatic.summarize;

import freemarker.template.Configuration;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import freemarker.template.Version;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class Summarizer {
    private static final Logger logger = LogManager.getLogger();

    private final String truncatePath;

    /** Map of data that will be used for the FreeMark template. */
    protected final Map<String, Object> templateData= new HashMap<>();
    /** FreeMarker configuration object. */
    protected final Configuration cfg;


    public Summarizer(String trunc) {
        this.truncatePath=trunc;
        this.cfg = new Configuration(new Version("2.3.23"));
        cfg.setDefaultEncoding("UTF-8");
        ClassLoader classLoader = Summarizer.class.getClassLoader();
        cfg.setClassLoaderForTemplateLoading(classLoader,"");
        parseTruncateData();
    }

    private void parseTruncateData() {
        logger.trace("Parsing the truncation data at {}", truncatePath);
        Map<String,String> truncatemap = new HashMap<>();
        try (
            BufferedReader br = new BufferedReader(new FileReader(truncatePath))) {
            String line;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                String[] fields = line.split(":");
                if (fields.length==2) {
                    truncatemap.put(fields[0],fields[1]);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace("Putting a total of {} items into the truncate map", truncatemap.size());
        templateData.put("truncate",truncatemap);
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

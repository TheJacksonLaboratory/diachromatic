package org.jax.diachromatic.align;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.FileNotFoundException;
import java.net.URL;

import static org.junit.jupiter.api.Assertions.assertNotNull;

public class DigestMapTest {


    private static String digestFilePath;

    @BeforeAll
    static void init() throws FileNotFoundException {
        ClassLoader classLoader = DigestMapTest.class.getClassLoader();
        URL url = classLoader.getResource("data/digestmap/hg19_digestedGenome.small.txt");
        if (url==null) {
            throw new FileNotFoundException("Could not find data/digestmap/hg19_digestGenome.small.txt");
        }
        digestFilePath = url.getFile();
    }

    @Test
    void testOpenFile() throws DiachromaticException  {
        DigestMap dmap = new DigestMap(digestFilePath);
        assertNotNull(dmap);
    }

}

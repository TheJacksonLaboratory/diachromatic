package org.jax.diachromatic.digest;

import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class RestrictionEnzymeTest {
    private static RestrictionEnzyme hindIII;
    private static List<RestrictionEnzyme> reList;

    @BeforeClass
    public static void setup() throws Exception {
        hindIII = new RestrictionEnzyme("HindIII", "A^AGCTT");
        ClassLoader classLoader = RestrictionEnzymeTest.class.getClassLoader();
        String restrictionEnzymesPath = classLoader.getResource("data/enzymelist.tab").getFile();
        reList=RestrictionEnzyme.parseRestrictionEnzymesFromFile(restrictionEnzymesPath);
    }
    @Test
    public void testGetName() {
        String expected="HindIII";
        assertEquals(expected, hindIII.getName());
    }
    @Test
    public void testGetRestrictionSite() {
        String expected="A^AGCTT";
        assertEquals(expected, hindIII.getSite());
    }

    /*
     * $ grep -v '#' enzymelist.tab | wc -l
     * 15
     */
    @Test
    public void testGetAllEnzymesFromFile() {
        assertEquals(15,reList.size());
    }

    @Test
    public void testContainsHindIII() {
        // this tests where the overridden "equals" method is working as expected
        // One HindIII object was created as above, and onbe was created from the file and placed into reList.
        assertTrue(reList.contains(hindIII));
    }

}

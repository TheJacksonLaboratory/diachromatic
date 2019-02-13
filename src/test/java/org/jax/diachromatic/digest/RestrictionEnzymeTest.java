package org.jax.diachromatic.digest;



import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class RestrictionEnzymeTest {
    private static RestrictionEnzyme hindIII;
    private static List<RestrictionEnzyme> reList;


    @BeforeAll
    static void setup() {
        hindIII = new RestrictionEnzyme("HindIII", "A^AGCTT");
        reList=RestrictionEnzyme.parseRestrictionEnzymes();
    }
    @Test
    void testGetName() {
        String expected="HindIII";
        assertEquals(expected, hindIII.getName());
    }
    @Test
    void testGetRestrictionSite() {
        String expected="A^AGCTT";
        assertEquals(expected, hindIII.getSite());
    }

    /*
     * $ grep -v '#' enzymelist.tab | wc -l
     * 15
     */
    @Test
    void testGetAllEnzymesFromFile() {
        assertEquals(15,reList.size());
    }

    @Test
    void testContainsHindIII() {
        // this tests where the overridden "equals" method is working as expected
        // One HindIII object was created as above, and onbe was created from the file and placed into reList.
        assertTrue(reList.contains(hindIII));
    }

}

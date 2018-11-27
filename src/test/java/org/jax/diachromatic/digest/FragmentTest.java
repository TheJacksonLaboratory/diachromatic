package org.jax.diachromatic.digest;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;


public class FragmentTest {


    @DisplayName("trying new parameterizied test function ")
    @ParameterizedTest(name = "{index} => pos={0}, samepos={1}")
    @CsvSource({
            "1, 1",
            "42, 42",
            "1042, 1042"
    })
    void parametrizedTest(int pos, int samepos) {
        int enyzymeNumber=19;
        Fragment f=new Fragment(enyzymeNumber,pos);
        assertEquals(samepos,f.position);
    }

    /** Test that we sort the {@link org.jax.diachromatic.digest.Fragment} objects in
     * ascending order of their position.
     */
    @Test
    void testFragmentSort() {
        int enyzymeNumber=19;
        Fragment f=new Fragment(enyzymeNumber,22);
        Fragment f2=new Fragment(enyzymeNumber,2);
        Fragment f3=new Fragment(enyzymeNumber,100);
        Fragment f4=new Fragment(enyzymeNumber,10);
        List<Fragment> fraglist = new ArrayList<>();
        fraglist.add(f);
        fraglist.add(f2);
        fraglist.add(f3);
        fraglist.add(f4);
        Collections.sort(fraglist);
        assertEquals(f2,fraglist.get(0));
        assertEquals(f4,fraglist.get(1));
        assertEquals(f,fraglist.get(2));
        assertEquals(f3,fraglist.get(3));
    }




}


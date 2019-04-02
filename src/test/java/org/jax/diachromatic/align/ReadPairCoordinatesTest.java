package org.jax.diachromatic.align;

import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

public class ReadPairCoordinatesTest {

    @Test
    void testEqualsFunction() {
        ReadPairCoordinates rpc1 = new ReadPairCoordinates(42,37,1);
        ReadPairCoordinates rpc2 = new ReadPairCoordinates(42,37,1);
        assertEquals(rpc1,rpc2);
    }

    @Test
    void testEqualsFunctionForDistinctRPC() {
        ReadPairCoordinates rpc1 = new ReadPairCoordinates(42,37,1);
        ReadPairCoordinates rpc2 = new ReadPairCoordinates(42,36,1);
        assertEquals(rpc1,rpc2);
    }

    @Test
    void testSet() {
        ReadPairCoordinates rpc1 = new ReadPairCoordinates(42,37,1);
        ReadPairCoordinates rpc2 = new ReadPairCoordinates(42,37,1);
        Set<ReadPairCoordinates> myset = new HashSet<>();
        myset.add(rpc1);
        myset.add(rpc2);
        assertEquals(1,myset.size());
    }
}

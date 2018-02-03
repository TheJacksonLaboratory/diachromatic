package org.jax.diachromatic.map;

/**
 * These are enumeration constants for errors that can be found by the Q/C of Capture Hi-C readpairs. The
 * quality checks are mainly performed in {@link ReadPair} and {@link SAMPairer}, see those classes for
 * details on the Q/C process and the meaning of these constants.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-06)
 */
public enum ErrorCode {

    READ_1_UNMAPPED("read 1 unmapped"),
    READ_2_UNMAPPED("read 2 unmapped"),
    READPAIR_UNMAPPED("readpair unmapped"),
    READ_1_MULTIMAPPED("read 1 multimapped"),
    READ_2_MULTIMAPPED("read 2 multimapped"),
    READPAIR_MULTIMAPPED("readpair multimapped"),
    INSERT_TOO_LONG("insert too long"),
    INSERT_TOO_SHORT("insert too short"),
    CIRCULARIZED_READ("circularized read (self-ligation)"),
    SAME_DANGLING_END("same dangling end"),
    CONTIGUOUS("contiguous"),
    RELIGATION("religation"),
    SAME_INTERNAL("same internal"),
    COULD_NOT_ASSIGN_TO_DIGEST("could not assign read pair to digest");

    private final String name;


    ErrorCode(String n) { name=n;}


    public String getName(){ return this.name; }




}

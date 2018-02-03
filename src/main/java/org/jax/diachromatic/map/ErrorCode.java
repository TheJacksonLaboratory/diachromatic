package org.jax.diachromatic.map;

public enum ErrorCode {

    READ_1_UNMAPPED("read 1 unmapped"),
    READ_2_UNMAPPED("read 2 unmapped"),
    READPAIR_UNMAPPED("readpair unmapped"),
    READ_1_MULTIMAPPED("read 1 multimapped"),
    READ_2_MULTIMAPPED("read 2 multimapped"),
    READPAIR_MULTIMAPPED("readpair multimapped");

    private final String name;


    ErrorCode(String n) { name=n;}


    public String getName(){ return this.name; }




}

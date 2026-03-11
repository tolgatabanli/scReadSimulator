package org.gobiws26;

public enum Config {
    INSTANCE;

    public static Integer READ_LENGTH;
    public static Integer FR_LENGTH;
    public static Integer SD;
    public static Integer TAIL_LENGTH;
    public static Double MUTATION_RATE;
    public static String OUTPUT_DIR; // will contain fw.fastq, rw.fastq and read.mappinginfo
}
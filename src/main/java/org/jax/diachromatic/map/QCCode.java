package org.jax.diachromatic.map;

/**
 * These are enumeration constants for errors that can be found by the Q/C of Capture Hi-C readpairs. The
 * quality checks are mainly performed in {@link ReadPair} and {@link SAMPairer}, see those classes for
 * details on the Q/C process and the meaning of these constants.
 * <p>
 * The following explanation of these artefacts is closely based on the HiCUP article: Wingett S,
 * Ewels P, Furlan-Magaril M et al. HiCUP: pipeline for mapping and processing Hi-C data. F1000Research 2015, 4:1310
 * </p>
 * <p>
 *     <b>Contiguous sequence artefacts</b>. Re-ligation or incomplete digestion leads to the generation of invalid contiguous
 *     sequences.
 *     <ul>
 *         <li>{@link #RELIGATION} occurs when adjacent ligation fragments re-attach. Adjacent fragments have the
 *         same orientation and thus the reads have opposite orientation. In our code, we determine whether the fragments are
 *         adjacent because their fragment numbers differ by 1. </li>
 *         <li>{@link #CONTIGUOUS} occurs when multiple fragments re-ligate and form a contig. Here, paired reads will not map to
 *              adjacent genomic restriction fragments (and thus fragment numbers differ by more than 1) the two reads
 *              are on different fragments that are not direct neighbors. However,  the
 *              paired reads are are located within one expected fragment size. This indicated thatthey are contiguous
 *              sequences that were not properly digested. The code demands that the contig size is above the lower
 *              threshold and below the upper threshold.</li>
 *
 *     </ul>
 * </p>
 * <p><b>Artefacts where the sequenced read-pair maps to a single restriction fragment</b>
 * Circularization by self-ligation by ligating to themselves. Such events are identifiable because the read pairs that
 * map to the same genomic restriction fragment are orientated away from each other when aligned to the reference genome.
 * This is tested for by the function {@code selfLigation()} in {@link ReadPair}. The function {@code danglingEnd},
 * also in {@link ReadPair}, tests if one of the read ends is within a threshold ({@code DANGLING_THRESHOLD}=7nt) of the
 * end of the restriction digest.
 *  <ul>
 *       <li>{@link #SAME_DANGLING_END} occurs if a readpair is not circularized but one of the reads overlaps with the
 *       end of a restriction digest.</li>
 *       <li>{@link #SAME_INTERNAL} occurs if a readpair is not circularized and neither of its ends overlaps
 *       with the end of the restriction digest.</li>
 *
 *  </ul>
 *
 * </p>
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-06)
 */
public enum QCCode {

    READ_1_UNMAPPED("read 1 unmapped"),
    READ_2_UNMAPPED("read 2 unmapped"),
    READPAIR_UNMAPPED("readpair unmapped"),
    READ_1_MULTIMAPPED("read 1 multimapped"),
    READ_2_MULTIMAPPED("read 2 multimapped"),
    READPAIR_MULTIMAPPED("readpair multimapped"),
    INSERT_TOO_LONG("insert too long"),
    INSERT_TOO_SHORT("insert too short"),
    RELIGATION("Re-ligation of adjacent restriction fragments"),
    SAME_DANGLING_END("same dangling end"),
    CONTIGUOUS("contiguous"),
    CIRCULARIZATION_INTERNAL("circularization (self-ligation) internal"),
    CIRCULARIZATION_DANGLING("circularization (self-ligation) dangling end"),
    SAME_INTERNAL("same internal"),
    COULD_NOT_ASSIGN_TO_DIGEST("could not assign read pair to digest");

    private final String name;


    QCCode(String n) { name=n;}


    public String getName(){ return this.name; }




}

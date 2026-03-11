package org.gobiws26;

import org.apache.commons.rng.simple.JDKRandomWrapper;
import org.apache.commons.statistics.distribution.BinomialDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Transcript {
    private static final Random RANDOM_SAMPLER = new Random(2025);
    private static final char[] ANTI_A = new char[] {'C', 'G', 'T'};
    private static final char[] ANTI_C = new char[] {'A', 'G', 'T'};
    private static final char[] ANTI_G = new char[] {'A', 'C', 'T'};
    private static final char[] ANTI_T = new char[] {'A', 'C', 'G'};

    private String chr;
    private String geneID;
    private String transcriptID;
    private Boolean isPositiveStranded = null;
    private int utrLength = 0;
    int readcount;
    

    private List<Region> exons = new ArrayList<>();
    private int[] cumulativeIntronLengths;
    private int[] cumulativeExonLengths;
    int length;


    int start;
    int end;

    public Transcript(String geneID, String transcriptID, int readcount) {
        this.geneID = geneID;
        this.transcriptID = transcriptID;
        this.readcount = readcount;

    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public void addExon(Region exon) {
        exons.add(exon);
        length += exon.length();
    }

    public void setStrand(Boolean isPositive) {
        this.isPositiveStranded = isPositive;
    }

    public void addToUTRLength(int utrLen) {
        this.utrLength += utrLen;
    }

    public int getUTRLength() {
        return utrLength;
    }


    public void simulateReads(FastaRandomAccessReader fastaReader, BufferedWriter mappingWriter, BufferedWriter rwFastq, int[] readid) throws IOException {
        exons.sort(Comparator.comparingInt(Region::start));
        int[] intronLengths = new int[exons.size() - 1];
        for (int intronIndex = 0; intronIndex < intronLengths.length; intronIndex++) {
            intronLengths[intronIndex] = exons.get(intronIndex + 1).start() - exons.get(intronIndex).end(); // ends are exclusive
        }
        this.cumulativeIntronLengths = new int[exons.size()]; // first element is Zero
        System.arraycopy(intronLengths, 0, this.cumulativeIntronLengths, 1, intronLengths.length);
        for (int intronIndex = 1; intronIndex < cumulativeIntronLengths.length; intronIndex++) {
            cumulativeIntronLengths[intronIndex] += cumulativeIntronLengths[intronIndex - 1];
        }
        this.cumulativeExonLengths = new int[exons.size()];
        this.cumulativeExonLengths[0] = exons.getFirst().length();
        for (int exonIndex = 1; exonIndex < cumulativeExonLengths.length; exonIndex++) {
            this.cumulativeExonLengths[exonIndex] = cumulativeExonLengths[exonIndex - 1] + exons.get(exonIndex).length();
        }


        // read the sequence and create reverse complement if negative stranded
        // iterate through exons and make transcript sequence
        char[] genomicSequence = fastaReader.fastaOf(this.chr, this.start, this.end);
        char[] transcriptSequence = new char[this.length];
        int transcriptIndex = 0;
        for (Region exon : exons) {
            for (int exonIndex = exon.start(); exonIndex < exon.end(); exonIndex++) {
                transcriptSequence[transcriptIndex++] = genomicSequence[exonIndex - this.start];
            }
        }
        if (!isPositiveStranded) {
            transcriptSequence = reverseComplement(transcriptSequence);
        }

        outerloop:
        for (int i = 0; i < readcount; i++) {
            String readIdNumber = String.valueOf(readid[0]++);
            mappingWriter.append(readIdNumber).append("\t").append(this.chr).append("\t").append(this.geneID).append("\t").append(this.transcriptID).append("\t");
            rwFastq.append("@").append(readIdNumber).append("\n");

            // Choose read start uniformly within the last bases of the transcript.
            // If no tail length provided, use the whole transcript. Otherwise,
            // sum three_prime_utr length and provided tail length to determine
            // how many bases from the transcript end are eligible.
            int tailLen;
            if (Config.TAIL_LENGTH == null) {
                tailLen = this.length;
            } else {
                tailLen = this.getUTRLength() + Config.TAIL_LENGTH;
                tailLen = Math.max(0, tailLen);
                tailLen = Math.min(this.length, tailLen);
            }
            int tailStart = Math.max(0, this.length - tailLen);

            // pick a start position inside tail such that read fits or is clipped at transcript end
            int maxPossibleStart = Math.max(tailStart, this.length - Config.READ_LENGTH);
            int readStart = tailStart;
            if (maxPossibleStart > tailStart) {
                readStart = tailStart + RANDOM_SAMPLER.nextInt(maxPossibleStart - tailStart + 1);
            }

            int readEnd = Math.min(readStart + Config.READ_LENGTH, this.length);
            Region readRegion = new Region(readStart, readEnd);

            // Record mapping for the single output read (`read2`)
            Region[] gReadVectors = getGenomicVectorsOf(readRegion);
            for (int j = 0; j < gReadVectors.length; j++) {
                mappingWriter.append(gReadVectors[j].toString());
                if (j < gReadVectors.length - 1) {
                    mappingWriter.append('|');
                }
            }
            mappingWriter.append('\t');
            mappingWriter.append(readRegion.toString()).append('\t');

            // read sequence (directly from transcript sequence; transcript was already oriented by strand)
            int readLen = readRegion.end() - readRegion.start();
            int numMutations = BinomialDistribution.of(readLen, Config.MUTATION_RATE / (double) 100)
                    .createSampler(new JDKRandomWrapper(new Random(2025))).sample();
            Set<Integer> mutPositionsSet = new HashSet<>();
            while (mutPositionsSet.size() < numMutations && readLen > 0) {
                mutPositionsSet.add(RANDOM_SAMPLER.nextInt(readLen));
            }
            int[] mutPositionsArr = mutPositionsSet.stream().mapToInt(num -> num).sorted().toArray();
            char[] mutSeq = introduceMutation(Arrays.copyOfRange(transcriptSequence, readRegion.start(), readRegion.end()), mutPositionsArr);
            // write read directly (no reverse-complement), because transcriptSEQ already reflects strand
            rwFastq.write(mutSeq);
            rwFastq.append("\n").append("+").append(readIdNumber).append("\n")
                    .append("I".repeat(Math.max(0, mutSeq.length))).append("\n");

            // mutation numbers
            for (int j = 0; j < mutPositionsArr.length; j++) {
                mappingWriter.append(String.valueOf(mutPositionsArr[j]));
                if (j < mutPositionsArr.length - 1) {
                    mappingWriter.append(',');
                }
            }
            mappingWriter.append("\n");
        }

    }

    // reads are in transcript coordinates, exons are in genomic.
    // Idea:    1) Convert read.start and read.end to genomic by using cumulative intron lengths:
    //              - Binary search with start/end
    //              - Use search index to get cumulative intron length -> amount of expansion
    //          2) split read into vectors by iterating exons
    private Region[] getGenomicVectorsOf(Region readInTranscript) {
        // if negative stranded, reverse the region
        if (!isPositiveStranded) {
            readInTranscript = new Region(this.length - readInTranscript.end(), this.length - readInTranscript.start());
        }

        // === 1 ===
        // how many introns are to the left
        int startIndex = Arrays.binarySearch(cumulativeExonLengths, readInTranscript.start() + 1);
        if (startIndex < 0) startIndex = -startIndex - 1;
        // how many are to the right
        int endIndex = Arrays.binarySearch(cumulativeExonLengths, readInTranscript.end());
        if (endIndex < 0) endIndex = -endIndex - 1;

        int numOfVectors = endIndex - startIndex + 1;

        // first elem of cumulativeIntronLengths is always zero
        int start = this.start + readInTranscript.start() + cumulativeIntronLengths[startIndex];
        int end = this.start + readInTranscript.end() + cumulativeIntronLengths[endIndex];

        // === 2 ===
        Region[] vectors = new Region[numOfVectors];
        int vecIdx = 0;
        for (int exonIndex = startIndex; exonIndex <= endIndex; exonIndex++) {
            Region exon = exons.get(exonIndex); // O(1) since ArrayList

            // if first region, leave the read start as is, and same for the end
            vectors[vecIdx++] = new Region(exonIndex == startIndex ? start : exon.start(),
                                            exonIndex == endIndex ? end : exon.end());
        }
        return vectors;
    }



    private char[] reverseComplement(char[] seq) {
        int n = seq.length;
        char[] revComp = new char[n];

        for (int i = 0; i < n; i++) {
            char c = seq[n - 1 - i];
            char comp = switch (c) {
                case 'A' -> 'T';
                case 'C' -> 'G';
                case 'G' -> 'C';
                case 'T' -> 'A';
                default -> throw new IllegalArgumentException("Cannot reverse-complement base: " + c);
            };
            revComp[i] = comp;
        }
        return revComp;
    }


    private char[] introduceMutation(char[] original, int[] positions) {
        for (int pos : positions) {
            int randomBaseIndex = RANDOM_SAMPLER.nextInt(3);
            original[pos] = switch (original[pos]) {
                case 'A' -> ANTI_A[randomBaseIndex];
                case 'C' -> ANTI_C[randomBaseIndex];
                case 'G' -> ANTI_G[randomBaseIndex];
                case 'T' -> ANTI_T[randomBaseIndex];
                default -> throw new IllegalArgumentException("Unexpected base: " + original[pos]);
            };
        }
        return original;
    }

    @Override
    public String toString() {
        return geneID + "\t" + transcriptID + "\t" + start + "\t" + end + "\t" + readcount + "\n" +
                "\t" + exons.stream().sorted(Comparator.comparingInt(Region::start)).map(Region::toString).collect(Collectors.joining("|"));
    }

}

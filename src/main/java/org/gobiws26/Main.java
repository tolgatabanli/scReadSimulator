package org.gobiws26;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.lang.management.ManagementFactory;


public class Main {
    static final Pattern patternTranscriptID = Pattern.compile("transcript_id\\s+\"([^\"]+)\"");

    public static void main(String[] args) throws IOException {
        File readcounts = null;
        File fasta = null;
        File fidx = null;
        File gtf = null;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-length":
                    if (i + 1 < args.length) {
                        Config.READ_LENGTH = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: Please give an integer after -length!");
                        return;
                    }
                    break;

                case "-tailLength":
                    if (i + 1 < args.length) {
                        Config.TAIL_LENGTH = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: Please give an integer after -tailLength!");
                        return;
                    }
                    break;

                case "-readcounts":
                    if (i + 1 < args.length) {
                        readcounts = new File(args[++i]);
                    } else {
                        System.err.println("Error: Please specify a TSV file after -readcounts!");
                        return;
                    }
                    break;

                case "-mutationrate":
                    if (i + 1 < args.length) {
                        Config.MUTATION_RATE = Double.valueOf(args[++i]);
                    } else {
                        System.err.println("Error: Please specify a number between 0 and 100 after -mutationrate!");
                        return;
                    }
                    break;

                case "-fasta":
                    if (i + 1 < args.length) {
                        fasta = new File(args[++i]);
                    } else {
                        System.err.println("Error: Please specify a FASTA file after -fasta!");
                        return;
                    }
                    break;

                case "-fidx":
                    if (i + 1 < args.length) {
                        fidx = new File(args[++i]);
                    } else {
                        System.err.println("Error: Please specify a FASTA file after -fidx!");
                        return;
                    }
                    break;

                case "-gtf":
                    if (i + 1 < args.length) {
                        gtf = new File(args[++i]);
                    } else {
                        System.err.println("Error: Please specify a GTF file after -gtf!");
                        return;
                    }
                    break;

                case "-od":
                    if (i + 1 < args.length) {
                        Config.OUTPUT_DIR = args[++i];
                    } else {
                        System.err.println("Error: Please specify an output directory after -od!");
                        return;
                    }
                    break;

                default:
                    System.err.println("Unknown argument: " + args[i]);
                    return;
            }
        }

        // Create fastaReader with FASTA and FAI
        FastaRandomAccessReader fastaReader = new FastaRandomAccessReader(fasta, fidx);

        BufferedReader lineCounter = new BufferedReader(new FileReader(readcounts));
        int numLines = 0;
        while (lineCounter.readLine() != null) numLines++;
        lineCounter.close();

        // ===== Read Counts ===== //
        Map<String, Transcript> transcripts = new HashMap<>((int) (numLines / 0.75 + 1)); // TranscriptID -> (TranscriptID, count, GeneID)
        BufferedReader reader = new BufferedReader(new FileReader(readcounts));
        String line;
        reader.readLine(); // skip the header
        while ((line = reader.readLine()) != null) {
            String[] fields = line.split("\t"); // GeneID - TranscriptID - Count
            transcripts.put(fields[1], new Transcript(fields[0], fields[1], Integer.parseInt(fields[2])));
        }

        // ===== GTF ===== //
        reader = new BufferedReader(new FileReader(gtf));
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("#")) continue;

            String[] fields = line.split("\t");
            String chr = fields[0];
            String feature = fields[2];
            int start = Integer.parseInt(fields[3]);
            int end = Integer.parseInt(fields[4]);
            Boolean isPositive = fields[6].equals("+");

                // accept transcript, exon, and three_prime_utr features
                if (!feature.equals("transcript") && !feature.equals("exon")
                    && !feature.equalsIgnoreCase("three_prime_utr")
                    && !feature.equalsIgnoreCase("three_prime_UTR")) continue;

            // extract transcript ID and check if it's in readcounts
            Matcher matcher = patternTranscriptID.matcher(fields[8]);
            String transcriptId;
            if (matcher.find()) {
                transcriptId = matcher.group(1);
                if (!transcripts.containsKey(transcriptId)) continue;
            } else continue;

            // set coordinates or record UTR length
            Transcript transcript = transcripts.get(transcriptId);
            if (feature.equals("transcript")) {
                transcript.setStart(start);
                transcript.setEnd(end);
                transcript.setChr(chr);
                transcript.setStrand(isPositive);
            } else if (feature.equals("exon")) {
                transcript.addExon(new Region(start, end + 1));
            } else if (feature.equalsIgnoreCase("three_prime_utr") || feature.equalsIgnoreCase("three_prime_UTR")) {
                // GTF coordinates are inclusive
                int utrLen = end - start + 1;
                transcript.addToUTRLength(utrLen);
            }
        }

        try {
            Files.createDirectories(Path.of(Config.OUTPUT_DIR));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File mappingInfo = new File(Config.OUTPUT_DIR + File.separator + "read.mappinginfo");
        File read2Fastq = new File(Config.OUTPUT_DIR + File.separator + "read2.fastq");

        Iterator<Map.Entry<String, Transcript>> transcriptIterator = transcripts.entrySet().iterator();
        int[] indexHolder = {0};

        try (
            BufferedWriter mapWriter = new BufferedWriter(new FileWriter(mappingInfo));
            BufferedWriter read2Writer = new BufferedWriter(new FileWriter(read2Fastq))
        ) {
            mapWriter.write("readid\tchr\tgene\ttranscript\tread2_regvec\tt_read2_regvec\tread2_mut\n");

            while (transcriptIterator.hasNext()) {
                Map.Entry<String, Transcript> entry = transcriptIterator.next();
                Transcript t = entry.getValue();

                t.simulateReads(fastaReader, mapWriter, read2Writer, indexHolder);
                // free memory hopefully!!
                transcriptIterator.remove();
            }
        }
    }
}

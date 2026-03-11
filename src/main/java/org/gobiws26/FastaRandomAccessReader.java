package org.gobiws26;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

public class FastaRandomAccessReader {
    File fasta;
    Map<String, long[]> faIndex;
    long time = 0;

    private static final CharsetDecoder decoder = StandardCharsets.UTF_8.newDecoder();


    public FastaRandomAccessReader(File fasta, File fidx) throws IOException {
        this.fasta = fasta;
        readFAI(fidx); // sets faIndex
        //faIndex.entrySet().stream().forEach(e -> System.out.println(e.getKey() + " " + Arrays.toString(e.getValue())));

    }

    private void readFAI(File file) throws IOException {
        faIndex = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;

            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");

                String chrKey = fields[0];
                long[] data = new long[4];

                // Position information
                for (int i = 0; i < 4; i++) {
                    data[i] = Long.parseLong(fields[i + 1]);
                }

                faIndex.put(chrKey, data);
            }
        }
    }

    // returns the FASTA sequence of a region defined by a chromosome number and start/end positions in the GTF
    public char[] fastaOf(String chr, int start, int end) throws IOException {
        long startTime = System.nanoTime();
        // ===== index settings ===== //
        //long chrLength = faIndex.get(chr)[0];
        long chrStart = faIndex.get(chr)[1];
        int basesPerLine = (int) faIndex.get(chr)[2];
        int bytesPerLine = (int) faIndex.get(chr)[3];

        // ===== coordinates ===== //
        start = start - 1; // indexing is 0-based while GTF is 1-based
        int transcriptLength = end - start; // now is end exclusive so no need for -1
        int lineSkipsBefore = start / basesPerLine;   // gives how many line skips before transcript starts
        int lineSkipsAfter = (start % basesPerLine + transcriptLength) / basesPerLine; // how many line skips the transcript has in the fasta
        long startInFile = chrStart + (long) lineSkipsBefore * bytesPerLine + (start % basesPerLine);

        // ===== Random Access ===== //
        try (RandomAccessFile raf = new RandomAccessFile(fasta, "r")) {
            raf.seek(startInFile);
            byte[] data = new byte[transcriptLength + lineSkipsAfter];
            int bytesRead = raf.read(data);

            if (bytesRead != -1) {
                //System.out.println("Chromosome: " + chr + "\tchrStart: " + chrStart + "\tstartInFile: " + startInFile + "\tlineSkipsBefore: " + lineSkipsBefore);
                CharBuffer charBuffer = decoder.decode(ByteBuffer.wrap(data, 0, bytesRead));
                char[] source = charBuffer.array();
                int length = charBuffer.length();

                // remove newline chars
                int writeIndex = 0;
                int readIndex = 0;
                while (readIndex < length) {
                    char c = source[readIndex++];
                    if (c != '\n') {
                        source[writeIndex++] = c;
                    }
                }

                // final array of correct size
                if (writeIndex == length) {
                    this.time += System.nanoTime() - startTime;
                    return source;
                }

                char[] cleanedArray = new char[writeIndex];
                System.arraycopy(source, 0, cleanedArray, 0, writeIndex);

                this.time += System.nanoTime() - startTime;
                return cleanedArray;
            }


        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        this.time += System.nanoTime() - startTime;
        return new char[0];
    }
}

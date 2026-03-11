package org.gobiws26;

public record Region(int start, int end) {

    public int length() {
        return end - start; // ends are exclusive
    }

    @Override
    public String toString() {
        return start + "-" + end;
    }
}

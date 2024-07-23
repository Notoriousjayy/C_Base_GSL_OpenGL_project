//
// Created by jorda on 7/23/2024.
//

#include "graph_io.h"

/* External variable definitions */
long io_errors = 0;

/* Static variables */
static char buffer[81];
static char* cur_pos = buffer;
static FILE* cur_file = NULL;

/* Static arrays and constants */
static char icode[256];
static long checksum_prime = (1L << 30) - 83;

static long magic;
static long line_no;
static long final_magic;
static long tot_lines;
static char more_data;

/* Static string for character mapping */
static char* imap = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_^~&@,;.:?!%#$+-*/|\\<=>()[]{}`'\" \n";

/* Static file name buffer */
static char file_name[20];


/* Function implementations */

/**
 * @brief Fill the buffer with the next line from the current file.
 */
static void fill_buf(void) {
    register char* p;
    if (!fgets(buffer, sizeof(buffer), cur_file)) {
        io_errors |= FILE_ENDED_PREMATURELY;
        buffer[0] = more_data = 0;
    }
    for (p = buffer; *p; p++);
    if (p-- == buffer || *p != '\n') {
        io_errors |= MISSING_NEWLINE;
        p++;
    }
    while (--p >= buffer && *p == ' ');
    *++p = '\n';
    *++p = 0;
    cur_pos = buffer;
}

/**
 * @brief Set up the icode array for character mapping.
 */
static void icode_setup(void) {
    register long k;
    register char* p;
    for (k = 0; k < 256; k++) {
        icode[k] = UNEXPECTED_CHAR;
    }
    for (p = imap, k = 0; *p; p++, k++) {
        icode[*p] = k;
    }
}

char imap_chr(long d) {
    return d < 0 || d > strlen(imap) ? '\0' : imap[d];
}

long imap_ord(char c) {
    if (!icode['1']) {
        icode_setup();
    }
    return (c < 0 || c > 255) ? UNEXPECTED_CHAR : icode[c];
}

long new_checksum(char* s, long old_checksum) {
    register long a = old_checksum;
    register char* p;
    for (p = s; *p; p++) {
        a = (a + a + imap_ord(*p)) % checksum_prime;
    }
    return a;
}

void graph_newline(void) {
    if (++line_no > tot_lines) more_data = 0;
    if (more_data) {
        fill_buf();
        if (buffer[0] != '*') {
            magic = new_checksum(buffer, magic);
        }
    }
}

long graph_eof(void) {
    return !more_data;
}

char graph_char(void) {
    if (*cur_pos) return (*cur_pos++);
    return '\n';
}

void graph_backup(void) {
    if (cur_pos > buffer) cur_pos--;
}

long graph_digit(char d) {
    icode[0] = d;
    if (imap_ord(*cur_pos) < d) return icode[*cur_pos++];
    return -1;
}

unsigned long graph_number(char d) {
    register unsigned long a = 0;
    icode[0] = d;
    while (imap_ord(*cur_pos) < d) {
        a = a * d + icode[*cur_pos++];
    }
    return a;
}

char str_buf[STR_BUF_LENGTH];

char* graph_string(char* p, char c) {
    while (*cur_pos && *cur_pos != c) {
        *p++ = *cur_pos++;
    }
    *p++ = 0;
    return p;
}

void graph_raw_open(char* f) {
    if (!icode['1']) {
        icode_setup();
    }
    cur_file = fopen(f, "r");

#ifdef DATA_DIRECTORY
    if (!cur_file && (strlen(DATA_DIRECTORY) + strlen(f) < STR_BUF_LENGTH)) {
        sprintf(str_buf, "%s%s", DATA_DIRECTORY, f);
        cur_file = fopen(str_buf, "r");
    }
#endif
    if (cur_file) {
        io_errors = 0;
        more_data = 1;
        line_no = magic = 0;
        tot_lines = 0x7fffffff;
        fill_buf();
    } else {
        io_errors = CANT_OPEN_FILE;
    }
}

long graph_open(char* f) {
    strncpy(file_name, f, sizeof(file_name) - 1);

    graph_raw_open(f);
    if (cur_file) {
        sprintf(str_buf, "* File \"%s\"", f);
        if (strncmp(buffer, str_buf, strlen(str_buf))) return (io_errors |= BAD_FIRST_LINE);

        fill_buf();
        if (*buffer != '*') return (io_errors |= BAD_SECOND_LINE);

        fill_buf();
        if (*buffer != '*') return (io_errors |= BAD_THIRD_LINE);

        fill_buf();
        if (strncmp(buffer, "* (Checksum parameters ", 23)) return (io_errors |= BAD_FOURTH_LINE);
        cur_pos += 23;
        tot_lines = graph_number(10);
        if (graph_char() != ',') return (io_errors |= BAD_FOURTH_LINE);
        final_magic = graph_number(10);
        if (graph_char() != ')') return (io_errors |= BAD_FOURTH_LINE);

        graph_newline();
    }
    return io_errors;
}

long graph_close(void) {
    if (!cur_file) return (io_errors |= NO_FILE_OPEN);
    fill_buf();
    sprintf(str_buf, "* End of file \"%s\"", file_name);
    if (strncmp(buffer, str_buf, strlen(str_buf))) io_errors |= BAD_LAST_LINE;
    more_data = buffer[0] = 0;

    if (fclose(cur_file) != 0) return (io_errors |= CANT_CLOSE_FILE);
    cur_file = NULL;
    if (line_no != tot_lines + 1) return (io_errors |= WRONG_NUMBER_OF_LINES);
    if (magic != final_magic) return (io_errors |= WRONG_CHECKSUM);
    return io_errors;
}

long graph_raw_close(void) {
    if (cur_file) {
        fclose(cur_file);
        more_data = buffer[0] = 0;
        cur_pos = buffer;
        cur_file = NULL;
    }
    return magic;
}

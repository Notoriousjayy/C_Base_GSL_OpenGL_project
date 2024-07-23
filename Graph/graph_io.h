#ifndef GB_IO_H
#define GB_IO_H

#include <stdio.h>
#include <string.h>

#ifdef SYSV
#include <string.h>
#else
#include <strings.h>
#endif

/* Constants */
#define UNEXPECTED_CHAR 127
#define STR_BUF_LENGTH 160

/* Error codes */
#define CANT_OPEN_FILE 0x1
#define CANT_CLOSE_FILE 0x2
#define BAD_FIRST_LINE 0x4
#define BAD_SECOND_LINE 0x8
#define BAD_THIRD_LINE 0x10
#define BAD_FOURTH_LINE 0x20
#define FILE_ENDED_PREMATURELY 0x40
#define MISSING_NEWLINE 0x80
#define WRONG_NUMBER_OF_LINES 0x100
#define WRONG_CHECKSUM 0x200
#define NO_FILE_OPEN 0x400
#define BAD_LAST_LINE 0x800

/* External variable declarations */
extern long io_errors;

/* Function declarations */
extern char imap_chr(long d);
extern long imap_ord(char c);
extern void graph_newline(void);
extern long new_checksum(char* s, long old_checksum);
extern long graph_eof(void);
extern char graph_char(void);
extern void graph_backup(void);
extern long graph_digit(char d);
extern unsigned long graph_number(char d);
extern char str_buf[STR_BUF_LENGTH];
extern char* graph_string(char* p, char c);
extern void graph_raw_open(char* f);
extern long graph_open(char* f);
extern long graph_close(void);

#endif /* GB_IO_H */

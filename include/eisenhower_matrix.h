//
// Created by jorda on 8/10/2024.
//

#ifndef EISENHOWER_MATRIX_H
#define EISENHOWER_MATRIX_H

typedef struct {
    char title[256];
    int priority;
    int quadrant;
} Task;

void eisenhower_load_tasks(const char *filename, Task **tasks, int *num_tasks);
void classify_tasks(Task *tasks, int num_tasks);
void plot_matrix(Task *tasks, int num_tasks);

#endif // EISENHOWER_MATRIX_H

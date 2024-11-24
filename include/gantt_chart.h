//
// Created by jorda on 8/10/2024.
//

#ifndef GANTT_CHART_H
#define GANTT_CHART_H

typedef struct Task {
    char title[256];
    int priority;
    int quadrant;
    char start_time[20];  // e.g., "2024-08-01"
    char end_time[20];    // e.g., "2024-08-10"
} Task;

void gantt_load_tasks(const char *filename, Task **tasks, int *num_tasks);
void plot_gantt_chart(Task *tasks, int num_tasks);

#endif // GANTT_CHART_H


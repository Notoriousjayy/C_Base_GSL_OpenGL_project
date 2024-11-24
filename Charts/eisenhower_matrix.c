//
// Created by jorda on 8/10/2024.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/eisenhower_matrix.h"

void eisenhower_load_tasks(const char *filename, Task **tasks, int *num_tasks) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[1024];
    int task_count = 0;
    while (fgets(line, sizeof(line), file)) {
        if (task_count == 0) { // Skip the header
            task_count++;
            continue;
        }

        Task task;
        sscanf(line, "%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%d,%*[^,],%[^\n]",
               &task.priority, task.title);

        (*tasks) = realloc(*tasks, sizeof(Task) * task_count);
        (*tasks)[task_count - 1] = task;
        task_count++;
    }

    *num_tasks = task_count - 1;
    fclose(file);
}

void classify_tasks(Task *tasks, int num_tasks) {
    for (int i = 0; i < num_tasks; i++) {
        if (tasks[i].priority == 1) {
            tasks[i].quadrant = 1;  // Quadrant 1: Urgent and Important
        } else if (tasks[i].priority == 2) {
            tasks[i].quadrant = 2;  // Quadrant 2: Not Urgent but Important
        } else if (tasks[i].priority == 3) {
            tasks[i].quadrant = 3;  // Quadrant 3: Urgent but Not Important
        } else {
            tasks[i].quadrant = 4;  // Quadrant 4: Neither Urgent nor Important
        }
    }
}

void plot_matrix(Task *tasks, int num_tasks) {
    FILE *gp = popen("gnuplot -persistent", "w");

    if (!gp) {
        fprintf(stderr, "Error opening GNPlot.\n");
        return;
    }

    // Setup GNPlot
    fprintf(gp, "set title 'Eisenhower Matrix'\n");
    fprintf(gp, "set xlabel 'Importance'\n");
    fprintf(gp, "set ylabel 'Urgency'\n");
    fprintf(gp, "set xrange [-1:2]\n");
    fprintf(gp, "set yrange [-1:2]\n");
    fprintf(gp, "set xtics ('Not Important' 0, 'Important' 1)\n");
    fprintf(gp, "set ytics ('Not Urgent' 0, 'Urgent' 1)\n");
    fprintf(gp, "set grid\n");

    // Plot the quadrants
    fprintf(gp, "plot '-' using 1:2:3 with labels offset 1,1 notitle\n");

    for (int i = 0; i < num_tasks; i++) {
        int x = 0, y = 0;
        switch (tasks[i].quadrant) {
            case 1:
                x = 1; y = 1; // Important and Urgent
                break;
            case 2:
                x = 1; y = 0; // Important but Not Urgent
                break;
            case 3:
                x = 0; y = 1; // Urgent but Not Important
                break;
            case 4:
                x = 0; y = 0; // Not Important and Not Urgent
                break;
        }

        fprintf(gp, "%d %d \"%s\"\n", x, y, tasks[i].title);
    }

    fprintf(gp, "e\n");
    fflush(gp);

    pclose(gp);
}

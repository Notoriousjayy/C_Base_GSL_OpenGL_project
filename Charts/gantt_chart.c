#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/gantt_chart.h"

void gantt_load_tasks(const char *filename, Task **tasks, int *num_tasks) {
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
        sscanf(line, "%*[^,],%[^,],%[^,],%[^,],%*[^,],%[^\n]",
               task.title, task.start_time, task.end_time);

        printf("Parsed Task: Title: %s, Start: %s, End: %s\n", task.title, task.start_time, task.end_time);

        *tasks = realloc(*tasks, sizeof(Task) * (task_count + 1));
        if (*tasks == NULL) {
            perror("Memory allocation error");
            exit(EXIT_FAILURE);
        }

        (*tasks)[task_count - 1] = task;
        task_count++;
    }

    *num_tasks = task_count - 1;
    fclose(file);
}


void plot_gantt_chart(Task *tasks, int num_tasks) {
    FILE *gp = popen("gnuplot -persistent", "w");

    if (!gp) {
        fprintf(stderr, "Error opening GNPlot.\n");
        return;
    }

    // Setup GNPlot
    fprintf(gp, "set title 'Gantt Chart'\n");
    fprintf(gp, "set xlabel 'Time'\n");
    fprintf(gp, "set ylabel 'Tasks'\n");
    fprintf(gp, "set xdata time\n");
    fprintf(gp, "set timefmt '%Y-%m-%d'\n");
    fprintf(gp, "set format x '%Y-%m-%d'\n");
    fprintf(gp, "set xtics rotate by -45\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set yrange [0:%d]\n", num_tasks + 1);
    fprintf(gp, "set ytics (");

    for (int i = 0; i < num_tasks; i++) {
        fprintf(gp, "\"%s\" %d", tasks[i].title, num_tasks - i);
        if (i < num_tasks - 1) {
            fprintf(gp, ", ");
        }
    }
    fprintf(gp, ")\n");

    // Plot the Gantt chart using horizontal bars
    fprintf(gp, "plot '-' using 2:1:3 with boxes title ''\n");

    for (int i = 0; i < num_tasks; i++) {
        fprintf(gp, "%d %s %s\n",
                num_tasks - i,
                tasks[i].start_time,
                tasks[i].end_time);
    }

    fprintf(gp, "e\n");
    fflush(gp);

    pclose(gp);
}

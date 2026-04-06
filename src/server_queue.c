/* 
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2024 Douglas Eleveld.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/// This file implements the server queue. 

/// It compiles to nothing when not compiling with the install server
/// option.
#include "buildflags.h"

#ifdef OPENPMX_SERVER
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <signal.h>

#include "defines.h"
#include "utils/c22.h"

#define QUEUE_DIR OPENPMX_INSTALL_PATH "/openpmx.queue"

static char global_ticket_path[PATH_MAX] = "";

static void cleanup(void)
{
	if (global_ticket_path[0]) {
		unlink(global_ticket_path);
		rmdir(QUEUE_DIR); /* only succeeds if we are the last */
	}
}

static int is_process_alive(pid_t pid)
{
	if (kill(pid, 0) == 0)
		return 1;
	return errno == EPERM;
}

static void make_ticket(char* ticket_name, size_t ticket_size)
{
	let my_pid = getpid();
	snprintf(ticket_name, ticket_size,
			 "%ld-%d", (long)time(NULL), my_pid);
	snprintf(global_ticket_path, sizeof(global_ticket_path),
			 "%s/%s", QUEUE_DIR, ticket_name);

	var f = fopen(global_ticket_path, "w");
	if (f) 
		fclose(f);
}

void server_queue(void)
{
	atexit(cleanup);

	/* make the queue directory */
	if (mkdir(QUEUE_DIR, 0755) && errno != EEXIST) {
		perror("mkdir");
		exit(EXIT_FAILURE);
	}

	/* Create unique ticket */
	char ticket_name[256];
	make_ticket(ticket_name, sizeof(ticket_name));

	let deadline = time(NULL) + OPENPMX_SERVER_TIMEOUT_S;
	int last_reported_rank = -1;
	int last_reported_tickets_found = -1;
	while (1) {
		if (time(NULL) >= deadline) {
			fprintf(stderr, "queue timed out after %d s\n",
					OPENPMX_SERVER_TIMEOUT_S);
			exit(EXIT_FAILURE);
		}

		/* get all the names of the running processes */
		struct dirent **namelist;
		int n = scandir(QUEUE_DIR, &namelist, NULL, alphasort);
		if (n < 0)
			break;

		/* look through entries looking for lowest rank and my rank */
		int my_rank = -1;
		char oldest_ticket[NAME_MAX + 1] = "";
		int tickets_found = 0;
		forcount(i, n) {
			if (namelist[i]->d_name[0] == '.')
				continue;

			if (tickets_found == 0)
				snprintf(oldest_ticket, sizeof(oldest_ticket),
						 "%s", namelist[i]->d_name);

			tickets_found++;

			if (strcmp(namelist[i]->d_name, ticket_name) == 0)
				my_rank = tickets_found;
		}
		/* cleanup namelist */
		forcount(i, n) 
			free(namelist[i]);
		free(namelist);

		/* remove stale leader */
		if (strcmp(oldest_ticket, ticket_name) != 0) {
			let hyphen = strrchr(oldest_ticket, '-');
			if (hyphen) {
				let leader_pid = atoi(hyphen + 1);
				if (!is_process_alive(leader_pid)) {
					char stale_path[PATH_MAX];
					snprintf(stale_path, sizeof(stale_path),
							 "%s/%s", QUEUE_DIR, oldest_ticket);
					unlink(stale_path);
					continue;
				}
			}
		}

		/* my ticket got lost, make and try again */	
		if (my_rank == -1) {
			make_ticket(ticket_name, sizeof(ticket_name));
			continue;
		}

		/* Its our turn now. The ticket should stay until we are done.
		 * If we exit and dont cleanup, a waiting process will remove
		 * our stale ticket */
		if (my_rank == 1)
			break;

		/* update the rank */
		if ((my_rank > 1 && my_rank != last_reported_rank) ||
			(last_reported_tickets_found != tickets_found)) {
			last_reported_rank = my_rank;
			last_reported_tickets_found = tickets_found;
		}
		
		sleep(2);
	}
}
#endif /* OPENPMX_SERVER */

#ifdef OPENPMX_NOTSERVER
void server_queue(void)
{
	return;
}
#endif /* OPENPMX_NOTSERVER */

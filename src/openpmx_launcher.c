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

/* Launcher: translate (if needed), compile, link and run an OpenPMX control
 * file. The following macros are supplied by the install script:
 * OPENPMX_CC, OPENPMX_CFLAGS, OPENPMX_LFLAGS, OPENPMX_LDLIBS, OPENPMX_PATH. */

#include "utils/c22.h"

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <limits.h>

#if !defined(OPENPMX_CC) || !defined(OPENPMX_CFLAGS) || !defined(OPENPMX_LFLAGS) || \
	!defined(OPENPMX_LDLIBS) || !defined(OPENPMX_PATH)
#error "OPENPMX_CC, OPENPMX_CFLAGS, OPENPMX_LFLAGS, OPENPMX_LDLIBS and OPENPMX_PATH must be defined"
#endif

#if defined(__MSYS__) || defined(__CYGWIN__) || defined(_WIN32)
#define EXEEXT ".exe"
#else
#define EXEEXT ".bin"
#endif

#if defined(_WIN32) && !defined(__MSYS__) && !defined(__CYGWIN__)
#define NATIVE_WIN32 1
#else
#define NATIVE_WIN32 0
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

static char *objfile = NULL;
static char *execfile = NULL;
static bool succeeded = false;

static void die_nomem(void)
{
	fprintf(stderr, "error: out of memory\n");
	exit(1);
}

/* printf into a freshly malloc'd string */
static char *xasprintf(const char *fmt, ...)
{
	va_list ap, ap2;
	va_start(ap, fmt);
	va_copy(ap2, ap);
	let n = vsnprintf(NULL, 0, fmt, ap);
	va_end(ap);
	if (n < 0)
		die_nomem();
	var s = mallocvar(char, n + 1);
	if (!s)
		die_nomem();
	vsnprintf(s, (size_t)n + 1, fmt, ap2);
	va_end(ap2);
	return s;
}

#if NATIVE_WIN32
/* quote a string for Windows cmd.exe using double quotes */
static char *shquote_win32(const char *s)
{
	let len = (int)strlen(s);
	var n = len * 2 + 3;
	var out = mallocvar(char, n);
	if (!out)
		die_nomem();
	var q = out;
	*q++ = '"';
	forcount(i, len) {
		if (s[i] == '"')
			*q++ = '\\';
		*q++ = s[i];
	}
	*q++ = '"';
	*q = 0;
	return out;
}
#else
/* quote a string for POSIX / MSYS2 / WSL / macOS shells with single quotes */
static char *shquote_posix(const char *s)
{
	let len = (int)strlen(s);
	var n = 2;
	forcount(i, len)
		n += (s[i] == '\'') ? 4 : 1;
	var out = mallocvar(char, n + 1);
	if (!out)
		die_nomem();
	var q = out;
	*q++ = '\'';
	forcount(i, len) {
		if (s[i] == '\'') {
			memcpyvar(q, "'\\''", 4);
			q += 4;
		} else
			*q++ = s[i];
	}
	*q++ = '\'';
	*q = 0;
	return out;
}
#endif

static char *shquote(const char *s)
{
#if NATIVE_WIN32
	return shquote_win32(s);
#else
	return shquote_posix(s);
#endif
}

/* the object and executable are removed, the .c file remains for
 * debugging purposes */
static void removefiles(void)
{
	static char **const files[] = { &objfile, &execfile };
	forarray(i, files)
		if (*files[i])
			remove(*files[i]);
}

static void cleanup(void)
{
	if (!succeeded) {
		printf("error: openpmx failed\n");
		fflush(stdout);
	}
	removefiles();
}

static void onsignal(int sig)
{
	(void)sig;
	removefiles();
	_exit(1);
}

/* run a command through the system shell, stderr merged into stdout */
static bool run(const char *cmd)
{
#if NATIVE_WIN32
	let full = xasprintf("(%s) 2>&1", cmd);
#else
	let full = xasprintf("( %s ) 2>&1", cmd);
#endif
	fflush(stdout);
	let rc = system(full);
	free(full);
	return rc == 0;
}

static bool has_path_sep(const char *path)
{
	return strchr(path, '/') != NULL || strchr(path, '\\') != NULL;
}

static const char *basename_of(const char *path)
{
	let s1 = strrchr(path, '/');
	let s2 = strrchr(path, '\\');
	var s = s1;
	if (!s || (s2 && s2 > s))
		s = s2;
	return s ? s + 1 : path;
}

int main(int argc, char **argv)
{
	if (argc < 2 || argc > 3) {
		fprintf(stderr, "usage: %s controlfile [datafile]\n", argc > 0 ? argv[0] : "openpmx");
		return 1;
	}
	let file_arg = argv[1];
	let data = (argc > 2) ? argv[2] : "";
	if (strlen(file_arg) >= PATH_MAX) {
		fprintf(stderr, "error: path too long\n");
		return 1;
	}

	/* control, object, and executable */
	objfile = xasprintf("%s.o", basename_of(file_arg));
	execfile = xasprintf("%s" EXEEXT, file_arg);

	atexit(cleanup);
	static const int sigs[] = { SIGINT, SIGTERM };
	forarray(i, sigs)
		signal(sigs[i], onsignal);

	var file = strdup(file_arg);
	if (!file)
		die_nomem();

	/* translate control file to C code */
	let dot = strrchr(basename_of(file), '.');
	if (!dot || !streq(dot + 1, "c")) {
		printf("openpmxtran %s %s\n", file, data);
		let qf = shquote(file);
		let qd = shquote(data);
		let cmd = xasprintf("%s/bin/openpmxtran %s %s", OPENPMX_PATH, qf, qd);
		let ok = run(cmd);
		free(cmd);
		free(qf);
		free(qd);
		if (!ok)
			return 1;
		let tmp = xasprintf("%s.c", file);
		free(file);
		file = tmp;
	}

	/* compile the C code, only show relative paths */
	char cwd[PATH_MAX];
	if (!getcwd(cwd, sizeof cwd))
		cwd[0] = 0;
	printf("compile %s\n", file);
	{
		let qf = shquote(file);
		let qo = shquote(objfile);
		let qprefix = shquote(cwd);
		let qinc = shquote(OPENPMX_PATH "/include");
		let cmd = xasprintf("%s %s -fdebug-prefix-map=%s=. -I%s -c %s -o %s",
		                    OPENPMX_CC, OPENPMX_CFLAGS, qprefix, qinc, qf, qo);
		let ok = run(cmd);
		free(cmd);
		free(qf);
		free(qo);
		free(qprefix);
		free(qinc);
		if (!ok)
			return 1;
	}

	/* we have object files, so link now */
	printf("link %s\n", objfile);
	{
		let qo = shquote(objfile);
		let qe = shquote(execfile);
		let qlib = shquote(OPENPMX_PATH "/lib");
		let cmd = xasprintf("%s %s -o %s -L%s %s -Wl,-Bstatic -lopenpmx "
		                    "-Wl,-Bdynamic -lgsl -lgslcblas -lm",
		                    OPENPMX_CC, qo, qe, qlib, OPENPMX_LFLAGS);
		let ok = run(cmd);
		free(cmd);
		free(qo);
		free(qe);
		free(qlib);
		if (!ok)
			return 1;
	}

	/* execute and clear exit message */
	printf("execute %s\n", file);
	{
#if NATIVE_WIN32
		let path = has_path_sep(execfile) ? strdup(execfile) : strdup(execfile);
#else
		let path = has_path_sep(execfile) ? strdup(execfile) : xasprintf("./%s", execfile);
#endif
		let qp = shquote(path);
		let ok = run(qp);
		free(qp);
		free(path);
		if (!ok)
			return 1;
	}

	succeeded = true;
	free(file);
	return 0;
}

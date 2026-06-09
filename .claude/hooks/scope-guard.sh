#!/usr/bin/env bash
#
# scope-guard.sh — PreToolUse guard for the SSC/ZFS project (branch ssc-zfs).
#
# Enforces the non-negotiable standing orders from CLAUDE.md mechanically, so a
# headless `claude -p` run cannot violate them by accident:
#
#   1. §8  benchmarks.md is the READ-ONLY numerical oracle. Never edit it to make
#          a test pass. -> hard DENY any write/edit.
#   2. §1  Z-vector / orbital-relaxed densities, response/relaxation terms, and
#          analytic gradients of D are OUT OF SCOPE. -> DENY edits to the source
#          files that implement them and tell the agent to log it in QUESTIONS.md
#          instead of silently expanding scope.
#
# Wired via .claude/settings.json on PreToolUse for Edit|Write|MultiEdit|NotebookEdit.
# Reads the tool-call JSON on stdin; emits a PreToolUse permission decision on stdout.
# Fails OPEN (allows) on any internal error so the guard can never wedge the build.

set -euo pipefail

# --- read the hook payload --------------------------------------------------
payload="$(cat || true)"

allow() { exit 0; }   # no output + exit 0 => defer to normal permission flow

deny() {
  # $1 = human-readable reason fed back to the model
  jq -n --arg r "$1" '{
    hookSpecificOutput: {
      hookEventName: "PreToolUse",
      permissionDecision: "deny",
      permissionDecisionReason: $r
    }
  }'
  exit 0
}

# If we cannot parse, do not block real work.
command -v jq >/dev/null 2>&1 || allow
[ -n "$payload" ] || allow

file_path="$(printf '%s' "$payload" | jq -r '.tool_input.file_path // .tool_input.notebook_path // empty' 2>/dev/null || true)"
[ -n "$file_path" ] || allow

base="$(basename "$file_path")"

# --- rule 1: benchmarks.md is read-only ground truth (CLAUDE.md §8) ----------
if [ "$base" = "benchmarks.md" ]; then
  deny "BLOCKED by scope-guard: benchmarks.md is the READ-ONLY numerical oracle (CLAUDE.md §8). \
Never edit reference values to make a test pass — fix the code, not the target. \
If a benchmark number is genuinely wrong, record it in QUESTIONS.md and stop."
fi

# --- rule 2: out-of-scope source files (CLAUDE.md §1, boundary in §5) --------
# Match on basename so the rule holds regardless of where the path is rooted.
case "$base" in
  tdhf_mrsf_z_vector.F90 | tdhf_mrsf_gradient.F90)
    deny "BLOCKED by scope-guard: $base implements Z-vector / orbital-relaxed densities / \
analytic-gradient machinery, which is explicitly OUT OF SCOPE for the first-order SS ZFS work \
(CLAUDE.md §1, boundary marker in §5). Do not pull these in. You may READ them for reference, \
but if the task seems to REQUIRE editing them, stop and record the scope conflict in QUESTIONS.md."
    ;;
esac

# --- default: allow ----------------------------------------------------------
allow

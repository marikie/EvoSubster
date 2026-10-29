# Implementation Agent Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Delegate approved implementation plans to a project-scoped GPT-6.1 Sol agent using medium reasoning effort.

**Architecture:** Define a standalone custom agent in `.codex/agents/implementer.toml`. Add parent orchestration instructions to the repository `AGENTS.md`: plan in the parent, delegate approved implementation, then review results. This is instruction-driven delegation, not a model-switch hook on the implementation popup.

**Tech Stack:** Codex custom-agent TOML and AGENTS.md.

**Spec:** User request in this conversation: configure an implementation subagent after planning.

## Global Constraints

- Scope configuration to this repository; preserve personal model settings and existing hooks/rules.
- Agent model: `gpt-6.1-sol`; `model_reasoning_effort = "medium"`.
- Inherit existing permissions; do not broaden sandbox or approval settings.
- Do not commit or push without a request.
- Respect repository instructions, existing user edits, and affected validation requirements.

## Review Focus

- Existing Codex binary may not support the current standalone agent format: check version and available configuration validation before claiming activation.
- Model availability: surface errors instead of silently selecting another model.
- Planning approval: do not implement while still in Plan mode or before the user approves the plan.
- Shared working files: only one implementation writer at a time; parallel editing requires separate worktrees.
- Configuration loading: verify agent discovery in a fresh session; distinguish static validation from a successful delegated run.

### Task 1: Configure and document the implementation agent

**Files:**
- Create: `.codex/agents/implementer.toml`
- Modify: `AGENTS.md`

**Interfaces:**
- Consumes: approved plan, scope, acceptance criteria, relevant test commands, working directory supplied by parent.
- Produces: implementation changes, executed validation results, limitations and deviations for parent review.

- [x] Check installed `codex --version` and configuration/help support against the official custom-agent documentation.
- [x] Create agent TOML with `name = "implementer"`, descriptive delegation guidance, model and effort above, and `developer_instructions` requiring scoped implementation, relevant tests, preservation of unrelated edits, and a concise handoff. Return material plan ambiguities to the parent; do not independently redesign the task or recursively delegate.
- [x] Add AGENTS.md instructions for the parent to delegate approved plans to `implementer`, supply complete context and test expectations, avoid simultaneous editing of its assigned files, and review its diff and validation before reporting completion. Keep the parent model unchanged. Do not promise automatic worktree isolation.
- [x] Parse the TOML with an available parser and inspect the resulting diff; preserve all existing hooks and configuration.
- [x] Reload in a fresh Codex session and confirm agent discovery. When execution is authorized, perform a read-only smoke delegation reporting the assigned role and working directory; verify model/effort from runtime metadata if available, not from agent self-report.
- [x] Report precisely which checks passed and whether activation/runtime verification remains pending. No repository pipeline is needed for configuration-only edits.

## Reference

Official schema and discovery locations: https://learn.chatgpt.com/docs/agent-configuration/subagents

## Execution Results

- Installed Codex: `0.162.0-alpha.17.2`.
- Added `.codex/agents/implementer.toml` and parent delegation instructions in `AGENTS.md`.
- TOML parsed successfully with the installed Python `toml` module; required fields and absence of permission overrides checked.
- Fresh Codex session discovered `implementer`. A persisted read-only CLI session spawned it with fresh context and completed the smoke test.
- Child session runtime `turn_context` independently confirmed `model = "gpt-6.1-sol"` and `effort = "medium"`; session metadata confirmed agent role `implementer`.
- An initial ephemeral smoke test failed because its parent rollout was unavailable; retrying in a persisted session succeeded.
- Independent read-only review: no important findings. `git diff --check`: passed.
- Existing hooks/rules, personal model settings, and unrelated changes were preserved. No commits or pushes.
- Ruling: edit configuration in the requested working directory rather than a separate worktree, so the active project can discover it; no application code was changed.
- `AGENTS.md` is locally ignored by Git in this checkout; its delegation rule is local unless deliberately shared later.

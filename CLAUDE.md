# Project instructions

## Markdown in this repo

Doug edits markdown here using RStudio's *Visual* mode, which rewrites the whole file to canonical form on every save. Anything written in a form Visual mode would reformat comes back on the next open as a diff nobody authored.

- **Do not hard-wrap.** Write each paragraph as a single line. Never add an `editor_options: markdown: wrap:` block to a file, and do not reinstate one that was removed.
- **No pipe tables.** Use bullet lists instead.
- **Keep inline code spans short.**
- **Escape markdown-significant characters in prose, the way Visual mode will.** Write `-\>` not `->`, and `\~` not `~`. Same for a literal `\<`, `\>` or `\|` in running text. Verified 2 August 2026: four unescaped lines were committed and RStudio rewrote all four on the next save, which is a phantom diff of exactly the kind these rules exist to prevent.
- **Do not escape inside code.** Backtick spans and fenced blocks are left verbatim by Visual mode, so `` `P = 2 -> 3` `` stays as it is. Escaping there would itself churn. Genuine markdown syntax is also untouched: `~~strikethrough~~` and a leading `>` blockquote marker stay as they are.
- **Re-read a file immediately before and after editing it**, since it may be open in the RStudio editor. If the result looks corrupted or duplicated, rewrite the whole file with Write rather than patching with Edit.
- **No em-dashes** anywhere, including commit messages. Use a plain double hyphen (`--`).

Wrapping "to N columns" is not a fix. It was tried and measured, and RStudio's canonicalisation cannot be reproduced outside RStudio, so no wrap width converges. Do not propose one.

`MarkdownWrap: None` in the `.Rproj` is what enforces this on the editor side.

The full rationale, the measurements behind it, and the rest of the project conventions are in `AGENTS.md`.

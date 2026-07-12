# Exporting the Agulhas Project from Claude to OpenAI Codex

A step-by-step guide for migrating the Agulhas project (files + conversation transcripts already downloaded from Claude.ai) into an OpenAI Codex project.

---

## Step 1: Organize your downloaded files into one folder

1. Create a folder on your computer, e.g. `Codex/Agulhas`
2. Drop in everything pulled from Claude:
   - Project knowledge files (docs, PDFs, notebooks, code)
   - Conversation transcripts / summaries
3. Optional but helpful — split into subfolders for clarity:
   ```
   Agulhas/
   ├── docs/
   ├── data/
   ├── code/
   └── transcripts/
   ```

---

## Step 2: Create an AGENTS.md file

This is Codex's equivalent of Claude's project custom instructions — it's read automatically at the start of every session in the folder.

1. In the root of `Codex/Agulhas`, create a plain text file named `AGENTS.md`
2. Paste in:
   - The project's goals and background (e.g. auxetic metamaterial research context, terminology, conventions)
   - Any standing instructions you gave Claude for this project
   - A pointer to the transcripts folder, e.g.:
     > See `/transcripts` for prior research discussion, decisions, and technical context from earlier work on this project.

---

## Step 3: Install Codex

Pick one:

- **Desktop app** — download from OpenAI's site, sign in with your ChatGPT account
- **CLI**:
  ```bash
  curl -fsSL https://chatgpt.com/codex/install.sh | sh      # macOS/Linux
  ```
  ```powershell
  irm https://chatgpt.com/codex/install.ps1 | iex           # Windows
  ```
- **VS Code extension** — install from the Extensions sidebar

---

## Step 4: Create the project in Codex

1. Open Codex (desktop app is easiest to start with)
2. Select **New thread** → create a **new project**
3. Point it at your `Codex/Agulhas` folder
4. Make sure **Local** ("Work locally") is selected so Codex works directly with your files on disk
5. Send a first check-in message:
   > "Read through the files in this project and summarize what you understand so far."
   This confirms everything loaded correctly before you dive into real tasks.

---

## Step 5: Feed in the transcript context

Transcripts are less structured than code, so give Codex an extra nudge:

- Reference the `/transcripts` folder explicitly in `AGENTS.md`, **or**
- Paste the most important summaries directly into your first message, so key decisions and context are in active context immediately rather than sitting in a file Codex might not read unprompted

---

## Note on working mode

Codex defaults to a code-and-repo workflow (running commands, editing files, opening PRs). If Agulhas is more research/writing-heavy than code-heavy:

- Use the desktop app in **local mode**
- You don't need to connect a GitHub repo unless you specifically want cloud/background task execution

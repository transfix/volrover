# heartbeat_job.py — the "hello world" of the volrover3 JobScheduler.
#
# Where the REPL runs one expression and returns, a *job* keeps running: the
# scheduler ticks its top-level `step(dt)` every frame (dt = seconds since the
# last tick), cooperatively, on the UI thread. That is the model for anything
# that unfolds over time — an animation, a simulation loop, a progress poll.
#
# HOW TO RUN (inside a running volrover3): open the "Python Console" dock, go to
# the "Jobs" tab, click "Load Script…", and pick this file. It appears in the
# jobs table and its Steps count climbs ~1/tick. Select it and click Stop to end
# it (or Interrupt to raise KeyboardInterrupt in it). See docs/EMBEDDED_PYTHON.md §12.
#
# THE CONTRACT: define a module-level `step(dt)`. Import-time code (below) runs
# ONCE when the job is submitted; `step` runs every tick. Raising from `step`
# stops the job (its status goes to error with the message shown in lastError),
# so a job "finishes" by raising — here we never do; the user stops it.

# ── import-time state (runs once, at submit) ─────────────────────────────────
_elapsed = 0.0  # total seconds this job has been ticking
_beats = 0      # how many one-second heartbeats we've emitted

# Optionally drive the live app: grab the real QMainWindow if we're running
# inside the GUI (headless job runs — e.g. tests — won't have one, so guard it).
try:
    import vrhost

    _window = vrhost.main_window()
    _title0 = _window.windowTitle()
except Exception:
    _window = None
    _title0 = ""

print("heartbeat_job: started — watch the Steps column climb; click Stop to end it.")


# ── per-tick work (runs every scheduler tick) ────────────────────────────────
def step(dt):
    global _elapsed, _beats
    _elapsed += dt
    # Emit a heartbeat once per whole second of accumulated time.
    if _elapsed >= _beats + 1:
        _beats += 1
        print(f"heartbeat_job: {_beats}s elapsed ({_elapsed:.2f}s real)")
        # If we have the live window, show the job driving the app over time by
        # tacking the beat count onto the title bar. This is the same live-C++
        # QMainWindow the menu_messagebox.py example mutates.
        if _window is not None:
            _window.setWindowTitle(f"{_title0}  —  heartbeat {_beats}s")

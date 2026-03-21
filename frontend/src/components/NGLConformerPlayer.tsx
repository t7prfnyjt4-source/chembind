import { useState, useEffect, useRef, useCallback } from "react";
import * as NGL from "ngl";

interface Atom {
  symbol: string;
  x: number;
  y: number;
  z: number;
}

interface Conformer {
  conf_id: number;
  energy: number;
  atoms: Atom[];
}

function conformerToSDF(atoms: Atom[]): string {
  const n = atoms.length;
  const lines = [
    "", "ChemBind Conformer", "",
    `${String(n).padStart(3)}${String(0).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`,
  ];
  for (const a of atoms) {
    lines.push(
      `${a.x.toFixed(4).padStart(10)}${a.y.toFixed(4).padStart(10)}${a.z.toFixed(4).padStart(10)} ${a.symbol.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`
    );
  }
  lines.push("M  END", "$$$$", "");
  return lines.join("\n");
}

export default function NGLConformerPlayer({
  conformers,
}: {
  conformers: Conformer[];
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);
  const [frame, setFrame] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [fps, setFps] = useState(5);
  const intervalRef = useRef<number | null>(null);

  // Initialize stage with first conformer
  useEffect(() => {
    if (!containerRef.current || !conformers.length) return;

    const stage = new NGL.Stage(containerRef.current, {
      backgroundColor: "white",
    });
    stageRef.current = stage;

    const sdf = conformerToSDF(conformers[0].atoms);
    const blob = new Blob([sdf], { type: "text/plain" });
    stage.loadFile(blob, { ext: "sdf" }).then((comp: any) => {
      comp.addRepresentation("ball+stick");
      comp.autoView();
    });

    return () => {
      stage.dispose();
      stageRef.current = null;
    };
  }, [conformers]);

  // Update on frame change
  useEffect(() => {
    if (!stageRef.current || !conformers.length) return;

    const stage = stageRef.current;
    stage.removeAllComponents();

    const sdf = conformerToSDF(conformers[frame].atoms);
    const blob = new Blob([sdf], { type: "text/plain" });
    stage.loadFile(blob, { ext: "sdf" }).then((comp: any) => {
      comp.addRepresentation("ball+stick");
      comp.autoView();
    });
  }, [frame, conformers]);

  const stopPlay = useCallback(() => {
    if (intervalRef.current !== null) {
      clearInterval(intervalRef.current);
      intervalRef.current = null;
    }
    setPlaying(false);
  }, []);

  useEffect(() => {
    if (!playing || !conformers.length) return;
    intervalRef.current = window.setInterval(() => {
      setFrame((f) => (f + 1) % conformers.length);
    }, 1000 / fps);
    return () => stopPlay();
  }, [playing, fps, conformers.length, stopPlay]);

  const current = conformers[frame];

  return (
    <div>
      <div style={{ fontWeight: 600, fontSize: 13, marginBottom: 4 }}>
        Conformer {frame + 1}/{conformers.length} — {current.energy.toFixed(2)} kcal/mol
      </div>

      <div
        ref={containerRef}
        style={{ width: 400, height: 400, border: "1px solid #555", borderRadius: 4 }}
      />

      <div style={{ display: "flex", gap: 8, alignItems: "center", marginTop: 8 }}>
        <button onClick={() => (playing ? stopPlay() : setPlaying(true))}>
          {playing ? "Pause" : "Play"}
        </button>
        <button onClick={() => { stopPlay(); setFrame(0); }}>Best</button>
        <input
          type="range"
          min={0}
          max={conformers.length - 1}
          value={frame}
          onChange={(e) => { stopPlay(); setFrame(Number(e.target.value)); }}
          style={{ flex: 1 }}
        />
      </div>

      <div style={{ display: "flex", gap: 8, alignItems: "center", marginTop: 4, fontSize: 13 }}>
        <label>FPS:</label>
        <input
          type="range" min={1} max={30} value={fps}
          onChange={(e) => setFps(Number(e.target.value))}
          style={{ width: 100 }}
        />
        <span>{fps}</span>
      </div>
    </div>
  );
}

import { useState, useEffect, useRef, useCallback } from "react";
import { apiFetch } from "../api/client";

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

interface ConformersResponse {
  conformers: Conformer[];
  cached: boolean;
}

const ELEMENT_COLORS: Record<string, string> = {
  C: "#333",
  N: "#3050F8",
  O: "#FF0D0D",
  S: "#FFFF30",
  P: "#FF8000",
  F: "#90E050",
  Cl: "#1FF01F",
  Br: "#A62929",
  I: "#940094",
};

function project2D(atoms: Atom[], width: number, height: number) {
  if (atoms.length === 0) return [];
  const xs = atoms.map((a) => a.x);
  const ys = atoms.map((a) => a.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const rangeX = maxX - minX || 1;
  const rangeY = maxY - minY || 1;
  const scale = Math.min((width - 40) / rangeX, (height - 40) / rangeY);
  const cx = width / 2, cy = height / 2;
  const midX = (minX + maxX) / 2, midY = (minY + maxY) / 2;

  return atoms.map((a) => ({
    sx: cx + (a.x - midX) * scale,
    sy: cy - (a.y - midY) * scale,
    symbol: a.symbol,
    color: ELEMENT_COLORS[a.symbol] || "#888",
  }));
}

export default function ConformerPlayer({ smiles }: { smiles: string }) {
  const [conformers, setConformers] = useState<Conformer[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [frame, setFrame] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [fps, setFps] = useState(5);
  const intervalRef = useRef<number | null>(null);

  const fetchConformers = async () => {
    setLoading(true);
    setError(null);
    try {
      const resp = await apiFetch<ConformersResponse>("/api/conformers", {
        method: "POST",
        body: { smiles, num_confs: 20 },
      });
      setConformers(resp.conformers);
      setFrame(0);
    } catch (e: any) {
      setError(e?.message || "Failed to generate conformers");
    } finally {
      setLoading(false);
    }
  };

  const stopPlay = useCallback(() => {
    if (intervalRef.current !== null) {
      clearInterval(intervalRef.current);
      intervalRef.current = null;
    }
    setPlaying(false);
  }, []);

  useEffect(() => {
    if (!playing || !conformers) return;
    intervalRef.current = window.setInterval(() => {
      setFrame((f) => (f + 1) % conformers.length);
    }, 1000 / fps);
    return () => stopPlay();
  }, [playing, fps, conformers, stopPlay]);

  if (!conformers) {
    return (
      <div>
        <button onClick={fetchConformers} disabled={loading}>
          {loading ? "Generating…" : "Generate Conformers"}
        </button>
        {error && <div style={{ color: "#c00", marginTop: 8, fontSize: 14 }}>{error}</div>}
      </div>
    );
  }

  const current = conformers[frame];
  const projected = project2D(current.atoms, 300, 300);

  return (
    <div style={{ marginTop: 12 }}>
      <div style={{ fontWeight: 600, marginBottom: 8 }}>
        Conformer {frame + 1}/{conformers.length} — Energy: {current.energy.toFixed(2)} kcal/mol
      </div>

      <svg width={300} height={300} style={{ border: "1px solid #555", borderRadius: 4, background: "#fafafa" }}>
        {projected.map((p, i) => (
          <g key={i}>
            <circle cx={p.sx} cy={p.sy} r={6} fill={p.color} />
            {p.symbol !== "C" && (
              <text x={p.sx} y={p.sy + 3} textAnchor="middle" fontSize={8} fill="#fff">{p.symbol}</text>
            )}
          </g>
        ))}
      </svg>

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
          type="range"
          min={1}
          max={30}
          value={fps}
          onChange={(e) => setFps(Number(e.target.value))}
          style={{ width: 100 }}
        />
        <span>{fps}</span>
      </div>
    </div>
  );
}

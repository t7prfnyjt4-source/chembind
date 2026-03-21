import { useState } from "react";
import { auth } from "../firebase/config";

const API_BASE = import.meta.env.VITE_API_BASE || "http://127.0.0.1:8000";

async function downloadExport(smiles: string, format: string) {
  const token = auth.currentUser ? await auth.currentUser.getIdToken() : null;
  const url = `${API_BASE}/api/export?smiles=${encodeURIComponent(smiles)}&format=${format}`;

  const res = await fetch(url, {
    headers: token ? { Authorization: `Bearer ${token}` } : {},
  });

  if (!res.ok) {
    const text = await res.text();
    throw new Error(`Export failed: ${res.status} ${text}`);
  }

  const blob = await res.blob();
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `molecule.${format}`;
  a.click();
  URL.revokeObjectURL(a.href);
}

export default function ExportButtons({ smiles }: { smiles: string }) {
  const [error, setError] = useState<string | null>(null);

  const handleExport = async (format: string) => {
    setError(null);
    try {
      await downloadExport(smiles, format);
    } catch (e: any) {
      setError(e?.message || "Export failed");
    }
  };

  return (
    <div style={{ marginTop: 8 }}>
      <div style={{ display: "flex", gap: 6 }}>
        <button onClick={() => handleExport("mol")} style={{ fontSize: 12 }}>Export MOL</button>
        <button onClick={() => handleExport("sdf")} style={{ fontSize: 12 }}>Export SDF</button>
        <button onClick={() => handleExport("cdxml")} style={{ fontSize: 12 }}>Export CDXML</button>
      </div>
      {error && <div style={{ color: "#c00", fontSize: 12, marginTop: 4 }}>{error}</div>}
    </div>
  );
}

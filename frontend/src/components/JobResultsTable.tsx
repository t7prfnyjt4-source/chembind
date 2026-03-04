import React from "react";
import type { JobItemDoc } from "../hooks/useJobItems";

function Badge({ ok }: { ok: boolean }) {
  return (
    <span
      style={{
        padding: "2px 8px",
        borderRadius: 999,
        background: ok ? "#ecfdf5" : "#fef2f2",
        color: ok ? "#065f46" : "#991b1b",
        fontWeight: 600,
        fontSize: 12,
      }}
    >
      {ok ? "Valid" : "Invalid"}
    </span>
  );
}

function fmt(n: unknown) {
  if (typeof n !== "number" || Number.isNaN(n)) return "—";
  return n.toFixed(3);
}

export function JobResultsTable({ items }: { items: Array<JobItemDoc & { id: string }> }) {
  return (
    <div style={{ marginTop: 12, overflowX: "auto", border: "1px solid #e5e7eb", borderRadius: 12 }}>
      <table style={{ width: "100%", borderCollapse: "collapse", minWidth: 900 }}>
        <thead>
          <tr style={{ background: "#f9fafb", textAlign: "left" }}>
            {["rowId", "SMILES", "Valid", "MW", "LogP", "HBD", "HBA", "tPSA", "Error"].map((h) => (
              <th
                key={h}
                style={{ padding: 10, borderBottom: "1px solid #e5e7eb", fontSize: 13, color: "#374151" }}
              >
                {h}
              </th>
            ))}
          </tr>
        </thead>

        <tbody>
          {items.map((it) => {
            const ok = Boolean(it.valid) && !it.error;
            return (
              <tr key={it.id} style={{ borderBottom: "1px solid #f3f4f6" }}>
                <td style={{ padding: 10, fontVariantNumeric: "tabular-nums" }}>{it.rowId ?? "—"}</td>
                <td style={{ padding: 10 }}>
                  <code style={{ fontSize: 12 }}>{it.smiles ?? "—"}</code>
                </td>
                <td style={{ padding: 10 }}>
                  <Badge ok={ok} />
                </td>
                <td style={{ padding: 10 }}>{fmt(it.mw)}</td>
                <td style={{ padding: 10 }}>{fmt(it.logp)}</td>
                <td style={{ padding: 10 }}>{typeof it.hbd === "number" ? it.hbd : "—"}</td>
                <td style={{ padding: 10 }}>{typeof it.hba === "number" ? it.hba : "—"}</td>
                <td style={{ padding: 10 }}>{fmt(it.tpsa)}</td>
                <td style={{ padding: 10, color: it.error ? "#991b1b" : "#6b7280" }}>{it.error ?? "—"}</td>
              </tr>
            );
          })}
        </tbody>
      </table>

      {items.length === 0 && <div style={{ padding: 12, color: "#6b7280" }}>No items yet…</div>}
    </div>
  );
}

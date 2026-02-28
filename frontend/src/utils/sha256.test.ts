import { describe, it, expect } from "vitest";
import { sha256Hex } from "./sha256";

describe("sha256Hex", () => {
  it("returns 64 hex chars", async () => {
    const hex = await sha256Hex("hello");
    expect(hex).toMatch(/^[0-9a-f]{64}$/);
  });

  it("is deterministic", async () => {
    const a = await sha256Hex("same");
    const b = await sha256Hex("same");
    expect(a).toBe(b);
  });

  it("changes when input changes", async () => {
    const a = await sha256Hex("a");
    const b = await sha256Hex("b");
    expect(a).not.toBe(b);
  });
});

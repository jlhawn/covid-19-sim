// The coordinate system is assumed to have its origin at the top-left.
class Rectangle {
	constructor(minX, minY, maxX, maxY) {
		this.minX = minX;
		this.minY = minY;
		this.maxX = maxX;
		this.maxY = maxY;
	}

	contains(point) {
		return this.minX <= point.x && point.x < this.maxX && this.minY <= point.y && point.y < this.maxY;
	}

	overlaps(rect) {
		// Simpler to determine if the rectangles *DO NOT* overlap and return logical-NOT of that.
		return !(this.maxX <= rect.minX || rect.maxX <= this.minX || this.maxY <= rect.minY || rect.maxY <= this.minY);
	}

	path2D() {
		var path = new Path2D();
		path.rect(this.minX, this.minY, this.maxX-this.minX, this.maxY-this.minY);
		return path;
	}
}

class EmptyValueList {
	get length() { return 0; }

	extend(next) { return next; }

	[Symbol.iterator]() {
		return {
			next() {
				return {done: true};
			}
		};
	}
}

const nilValueList = new EmptyValueList();

class ValueList {
	constructor(value) {
		this.head = {
			value: value,
			next: nilValueList
		};
		this.tail = this.head;
		this.size = 1;
	}

	get length() { return this.size; }

	extend(otherList) {
		if (otherList !== nilValueList) {
			this.tail.next = otherList.head;
			this.tail = otherList.tail;
			this.size += otherList.size;
		}
		return this;
	}

	[Symbol.iterator]() {
		return {
			ref: this.head,
			item: {value: null, done: false},
			next() {
				if (this.ref === nilValueList) {
					return {done: true};
				}

				this.item.value = this.ref.value;
				this.ref = this.ref.next;
				return this.item;
			}
		};
	}
}

class QuadTreeItem {
	constructor(point, ref, node) {
		this.point = point;
		this.ref = ref;
		this.node = node;
	}

	move(point) {
		if (this.node.contains(point)) {
			return; // Still in this node.
		}

		var newItem = null;
		var ancestor = this.node.parent;
		while (ancestor !== null) {
			if (ancestor.contains(point)) {
				// Insert and get new item;
				newItem = ancestor.insert(point, this.ref);
				break;
			}
			ancestor = ancestor.parent;
		}

		// Remove this item from its current node.
		this.remove();
		if (newItem === null) {
			console.log("No newItem", this);
			throw new Error("no new item");
		}
		// Attach this item to its new node.
		this.node = newItem.node;
		this.node.item = this;
	}

	remove() {
		this.node.item = null;

		var ancestor = this.node.parent;
		while (ancestor !== null && ancestor.simplify()) {
			ancestor = ancestor.parent;
		}
	}
}

class QuadTree extends Rectangle {
	constructor(minX, minY, maxX, maxY, parent, quadrant) {
		super(minX, minY, maxX, maxY);

		this.parent = parent || null;
		this.id = (parent && parent.id || "") + (quadrant || "");
		this.item = null;
		this.subTrees = null;
	}

	isEmpty() {
		return this.item === null && this.subTrees === null;
	}

	simplify() {
		if (this.subTrees === null) {
			return false; // Nothing to simplify.
		}

		// If any of the subtrees have any other subtrees then
		// we can't simplify either. That implies that there are
		// multiple items in that subtree.
		var directSubItemCount = 0;
		var directSubItem = null;
		for (var subTree of this.subTrees) {
			if (subTree.subTrees !== null) {
				return false; // subtree with subtrees.
			}
			if (subTree.item !== null) {
				directSubItemCount++;
				directSubItem = subTree.item;
			}
		}

		if (directSubItemCount === 0) {
			// This is a bug.
			throw new Error("all subtrees have no item or sub-items: "+this.id);
		}

		if (directSubItemCount !== 1) {
			return false; // Can only simplify if there is a single sub-item;
		}

		this.item = directSubItem;
		directSubItem.node = this;
		this.subTrees = null;
		return true;
	}	

	split() {
		var midX = (this.minX + this.maxX) / 2.0;
		var midY = (this.minY + this.maxY) / 2.0;

		this.subTrees = [
			new QuadTree(this.minX, this.minY, midX, midY, this, "A"), // Top-Left.
			new QuadTree(midX, this.minY, this.maxX, midY, this, "B"), // Top-Right.
			new QuadTree(this.minX, midY, midX, this.maxY, this, "C"), // Bottom-Left.
			new QuadTree(midX, midY, this.maxX, this.maxY, this, "D"), // Bottom-Right.
		];

		for(var subTree of this.subTrees) {
			if (subTree.contains(this.item.point)) {
				subTree.item = this.item;
				subTree.item.node = subTree;
				this.item = null;
				return;
			}
		}
	}

	insert(point, ref) {
		if (this.isEmpty()) {
			this.item = new QuadTreeItem(point, ref, this);
			return this.item;
		}

		if (this.subTrees === null) {
			this.split();
		}

		for (var subTree of this.subTrees) {
			if (subTree.contains(point)) {
				return subTree.insert(point, ref);
			}
		}
	}

	itemsInRange(rect) {
		if (this.isEmpty() || !this.overlaps(rect)) {
			return nilValueList;
		}

		if (this.subTrees === null) {
			if (rect.contains(this.item.point)) {
				return new ValueList(this.item.ref);
			}
			return nilValueList;
		}

		var vals = nilValueList;
		for (var subTree of this.subTrees) {
			vals = vals.extend(subTree.itemsInRange(rect));
		}
		return vals;
	}

	path2Ds() {
		var paths = [this.path2D()];
		if (this.subTrees !== null) {
			for (var subTree of this.subTrees) {
				paths.push(...subTree.path2Ds());
			}
		}
		return paths;
	}
}

class Point {
	constructor (x, y) {
		this.x = x;
		this.y = y;
	}

	translate(vec) {
		this.x += vec.dx;
		this.y += vec.dy;
	}

	vectorFrom(point) {
		return new Vector(this.x - point.x, this.y - point.y);
	}
}

class Circle extends Point {
	constructor(x, y, r) {
		super(x, y);

		this.r = r || 10;
	}

	boundingRect(padding) {
		var padding = padding || 0;
		var extent = this.r + padding;
		return new Rectangle(this.x - extent, this.y - extent, this.x + extent, this.y + extent);
	}

	overlaps(circle) {
		var dx = this.x - circle.x;
		var dy = this.y - circle.y;
		var rR = this.r + circle.r;
		var xSqrd = dx * dx;
		var ySqrd = dy * dy;
		var rSqrd = rR * rR;
		return (xSqrd + ySqrd) < rSqrd;
	}

	hitsWall(bWall) {
		// Really want to avoid doing sqrt.
		var rSqrd = this.r * this.r;

		var num = (bWall.a * this.x) + (bWall.b * this.y) + bWall.c;
		num *= num;

		var den = bWall.a * bWall.a + bWall.b * bWall.b;

		return (num/den) < rSqrd;
	}

	path2D() {
		var path = new Path2D();
		path.arc(this.x, this.y, this.r, 0, 2 * Math.PI);
		return path;
	}
}

class Vector {
	constructor(dx, dy) {
		this.dx = dx;
		this.dy = dy;
	}

	scale(factor) {
		return new Vector(this.dx *factor, this.dy * factor);
	}

	scalarProduct(vec) {
		return this.dx * vec.dx  + this.dy * vec.dy;
	}

	magnitude() {
		return Math.sqrt(this.scalarProduct(this));
	}

	add(vec) {
		return new Vector(this.dx + vec.dx, this.dy + vec.dy);
	}

	subtract(vec) {
		return new Vector(this.dx - vec.dx, this.dy - vec.dy);
	}

	projectOnto(vec) {
		// Let Bu be a unit vector in the direction of B, then
		// A projected onto B = A.B/|B| * Bu
		// Bu is equal to B/|B|   so, the projection becomes
		// A.B/|B| * B/|B| -> A.B/(|B|*|B|) * B
		// note that B.B = |B|*|B|
		// So the projection simplifies to (A.B)/(B.B)*B
		return vec.scale(this.scalarProduct(vec)/vec.scalarProduct(vec));
	}

	reflect(nVec) {
		// Reflect this vector against a surface which is defined by the
		// given vector which is normal to that surface. The normal vector
		// need not be a unit vector.
		return this.subtract(this.projectOnto(nVec).scale(2))
	}

	unitVector() {
		return this.scale(1/this.magnitude());
	}
}

function inelasticCollision(x1, x2, v1, v2, m1, m2, c) {
	// Normal vector from x2 to x1. (not necessarily a unit vector!)
	var nVec = x1.vectorFrom(x2);

	var v1Proj = v1.projectOnto(nVec);
	var v2Proj = v2.projectOnto(nVec);

	var vDiff = v2Proj.subtract(v1Proj);

	var vNet = vDiff.scale(c*m2).add(v1Proj.scale(m1).add(v2Proj.scale(m2))).scale(1/(m1+m2));

	return v1.subtract(v1Proj).add(vNet);
}

function elasticCollision(x1, x2, v1, v2, m1, m2) {
	// Normal vector from x2 to x1. (not necessarily a unit vector!)
	var nVec = x1.vectorFrom(x2);
	// Velocity difference vector.
	var vVec = v1.subtract(v2);

	return v1.subtract(vVec.projectOnto(nVec).scale(2*m2/(m1+m2)));
}

const healthy = 0;
const infected = 1;
const recovered = 2;
const infectionDuration = 200;

class BouncingBall {
	// A: x(t) = dx1 * t + x1
	//    y(t) = dy1 * t + y1
	//
	// B: x(t) = dx2 * t + x2
	//    y(t) = dy2 * t + y2
	//
	// B-A: x(t) = (dx2-dx1) * t + (x2-x1)
	//      y(t) = (dy2-dy1) * t + (y2-y1)

	constructor(qTree, circle, vector, density) {
		this.qTreeItem = qTree.insert(circle, this);
		this.circle = circle;
		this.vector = vector || (new Vector(0, 0));
		this.mass = (density || 1) * circle.r * circle.r; // Proportional to area.

		this.id = "b" + BouncingBall.nextID++;

		this.health = healthy;
		this.timeSinceInfection = null;
	}

	boundingRect(buffer) {
		return this.circle.boundingRect(buffer);
	}

	overlaps(bBall) {
		return this.circle.overlaps(bBall.circle);
	}

	hitsWall(wall) {
		return this.circle.hitsWall(wall);
	}

	move(dt) {
		this.circle.translate(this.vector.scale(dt));
		this.qTreeItem.move(this.circle);
		this.maybeRecover(dt);
	}

	maybeInfect(bBall) {
		if (this.health === infected && bBall.health < infected) {
			bBall.health = infected;
			bBall.timeSinceInfection = 0;
		}
	}

	maybeRecover(dt) {
		this.timeSinceInfection += dt;
		if (this.health === infected && this.timeSinceInfection > infectionDuration) {
			this.health = recovered;
		}
	}

	collideWith(bBall) {
		if (bBall instanceof BoundaryWall) {
			return bBall.collideWith(this);
		}

		var v1 = this.vector;
		var v2 = bBall.vector;

		var m1 = this.mass;
		var m2 = bBall.mass;

		var x1 = this.circle;
		var x2 = bBall.circle;

		this.vector = elasticCollision(x1, x2, v1, v2, m1, m2);
		bBall.vector = elasticCollision(x2, x1, v2, v1, m2, m1);

		this.maybeInfect(bBall);
		bBall.maybeInfect(this);
	}

	timeUntilCollision(other) {
		if (other instanceof BoundaryWall) {
			return this.timeUntilWallCollision(other);
		}
		return this.timeUntilBallCollision(other);
	}

	timeUntilBallCollision(bBall) {
		// We can quickly tell if the balls may collide based on whether the
		// current trajectory of one ball is toward the other from its own
		// reference frame (i.e., its speed is zero).
		// This is a normal vector from this ball to the other.
		var nVec = bBall.circle.vectorFrom(this.circle);
		// This is a vector for the relative velocity of this ball.
		var relVelVec = this.vector.subtract(bBall.vector);
		// The magnitude of the relative velocity vector projected onto the
		// normal vector is proportional to scalar product of those two vectors.
		// If the scalar product is negative then this ball is on a trajectory
		// away from the other ball.
		if (relVelVec.scalarProduct(nVec) < 0) {
			// console.log("heuristic is true");
			return null;
		}

		var dx = this.vector.dx - bBall.vector.dx;
		var dy = this.vector.dy - bBall.vector.dy;

		var x = this.circle.x - bBall.circle.x;
		var y = this.circle.y - bBall.circle.y;

		var d = this.circle.r + bBall.circle.r;

		var a = dx*dx + dy*dy;
		var b = 2 * (dx*x + dy*y);
		var c = x*x + y*y - d*d;

		if (c < 0) {
			// console.log("DEFORMATION: c="+c);
			// var overlap = d - Math.sqrt(x*x + y*y);

			// var d1 = this.mass/(this.circle.r*this.circle.r);
			// var d2 = bBall.mass/(bBall.circle.r*bBall.circle.r);

			// this.circle.r -= 2*overlap*(d2/(d1+d2));
			// bBall.circle.r -= 2*overlap*(d1/(d1+d2));

			// return null;

			console.log({dx, dy, x, y, a, b, c, d});
			throw new Error("c is "+c+" implying that "+this.id+" and "+bBall.id+" already overlap");
		}

		var discriminant = b*b - 4*a*c;

		// If the discriminant is zero, then the balls just slide by one
		// another (not transferring energy). If it is negative then the
		// balls do not collide at all (complex/imaginary result).
		if (discriminant <= 0) {
			return null
		}

		var dscrSqrt = Math.sqrt(discriminant);

		var t1 = (-b + dscrSqrt)/(2*a);
		var t2 = (-b - dscrSqrt)/(2*a);

		// console.log("t1: ", t1);
		// console.log("t2: ", t2);

		if (t1 < 0 || t2 < 0) {
			// Collision is in the past.
			return null;
		}

		return (t1 < t2) ? t1 : t2;
	}

	timeUntilWallCollision(bWall) {
		// Once we have the distance from the edge of the ball to the wall,
		// we need to find the magnitude of the the ball's velocity vector
		// projected onto the wall's normal unit vector. Because it's being
		// projected onto a unit vector, the magnitude of the projected vector
		// can be calculated as the absolute value of the scalar product. If
		// the projected vector is going away from the wall then the scalar
		// product will be positive, meaning we do not expect a collision.
		var relativeVelocity = this.vector.scalarProduct(bWall.nVec);
		if (relativeVelocity >= 0) {
			return null;
		}

		// The position vector from the wall's defining point to the center of
		// the ball.
		var pVec = this.circle.vectorFrom(bWall.point);
		// Projecting this onto the wall's normal vector gives us the position
		// vector from the point on the wall closest to the ball to the center
		// of the ball. The wall's normal vector is a unit vector so the
		// magnitude of this vector is the scalar product. Subtract the ball's
		// radius to get distance to the edge of the ball.
		var d = pVec.scalarProduct(bWall.nVec) - this.circle.r;
		if (d < 0) {
			console.log("Ball Over Wall", this, bWall);
			throw new Error("ball over wall");
		}

		// Dividing by magnidute of ball's velocity along this path gives the
		// time until the collision.
		return Math.abs(d/relativeVelocity);
	}

	vectorPath2D() {
		var path = new Path2D();
		path.moveTo(this.circle.x, this.circle.y);
		path.lineTo(this.circle.x + 2*this.vector.dx, this.circle.y + 2*this.vector.dy);
		return path;
	}

	draw(ctx) {
		var fillStyle = "white";
		switch (this.health) {
		case healthy:
			fillStyle = "skyblue";
			break;
		case infected:
			fillStyle = "orange";
			break;
		case recovered:
			fillStyle = "mediumpurple";
			break;
		}

		ctx.fillStyle = fillStyle;
		ctx.fill(this.circle.path2D())
	}
}

BouncingBall.nextID = 0;

class BoundaryWall {
	constructor(point, nVec, region) {
		// point(a, b), nVec<A, B>
		// Equation: A(x - a) + B(y - b) = 0
		// Simplifies to Ax + By + (-Aa - Bb) = 0
		this.nVec = nVec.unitVector();
		this.point = point;
		this.a = nVec.dx;
		this.b = nVec.dy;
		this.c = -this.a*point.x - this.b*point.y;

		this.region = region;

		this.id = "w" + BoundaryWall.nextID++;
	}

	collideWith(bBall) {
		bBall.vector = bBall.vector.reflect(this.nVec);
	}
}

BoundaryWall.nextID = 0;

class CollisionEvent {
	constructor(t, b1, b2) {
		this.t = t; // Time when this event occurs. NOT current time.
		if (b1.id < b2.id) {
			this.b1 = b1; // First Belligerent.
			this.b2 = b2; // Second Belligerent.
		} else {
			this.b1 = b2;
			this.b2 = b1;
		}

		this.key = [this.t, this.b1.id, this.b2.id]
		this._balls = [];
		this._ballIDs = [];
		for (var b of [this.b1, this.b2]) {
			if (b instanceof BouncingBall) {
				this._balls.push(b);
				this._ballIDs.push(b.id);
			}
		}

		this.isMaybe = false;
		this.index = -1;
	}

	objectsEqual(other) {
		return this.b1.id === other.b1.id && this.b2.id === other.b2.id;
	}

	balls() {
		return this._balls;
	}

	ballIDs() {
		return this._ballIDs;
	}

	perform() {
		this.b1.collideWith(this.b2);
	}

	reevaluate(t, frameDuration) {
		var dt = this.b1.timeUntilCollision(this.b2);
		if (dt === null) {
			return false; // Event will not happen.
		}
		if (t + dt > frameDuration) {
			return false; // Event occurs after this frame ends.
		}
		this.t = t + dt;
		return true;
	}

	lessThan(event) {
		for (var i = 0; i < 3; i++) {
			if (this.key[i] < event.key[i]) {
				return true;
			}
			if (this.key[i] > event.key[i]) {
				return false;
			}
		}
		return false;
	}

	lessThanOrEqual(event) {
		for (var i = 0; i < 3; i++) {
			if (this.key[i] < event.key[i]) {
				return true;
			}
			if (this.key[i] > event.key[i]) {
				return false;
			}
		}
		return true;
	}
}

function parentIndexOf(idx) {
	return Math.floor((idx+1)/2) - 1;
}

function leftChildIndexOf(idx) {
	return idx*2 + 1;
}

function rightChildIndexOf(idx) {
	return leftChildIndexOf(idx)+1;
}

class CollisionEventsMinHeap {
	constructor() {
		this.arr = [];

		// Maps a ball ID to a set of events which involve that ball.
		// The event object will have its index set as a field.
		this.ballEvents = new Map();
	}

	check(frameNum) {
		// Each item of the array should have a time which is
		// greater than or equal to the last. The first event
		// should not have a time less than zero.
		for (var idx = 0; idx < this.arr.length; idx++) {
			var parentIdx = parentIndexOf(idx);
			if (idx === 0 && this.arr[idx].t < 0 || idx > 0 && this.arr[parentIdx].t > this.arr[idx].t) {
				throw new Error("invalid heap property");
			}
		}
	}

	isEmpty() {
		return this.arr.length === 0;
	}

	add(event) {
		// Check for other events with the same balls.
		// There may be an existing event for either ball involved.
		// If there is an existing event which involves both balls,
		// keep only the earliest one. If there is an event involving
		// only one ball, keep the earlier one but mark the later one
		// as "maybe" so that it can be re-evaluated later. The earlier
		// event will affect the outcome of the later event.
		var existingEvents = new Set();
		for (var id of event.ballIDs()) {
			var ballEvents = this.ballEvents.get(id);
			if (!ballEvents) { continue; } // No existing event for this ball ID.
			for (var ballEvent of ballEvents) {
				existingEvents.add(ballEvent);
			}
		}

		for (var existingEvent of existingEvents) {
			if (existingEvent.objectsEqual(event)) {
				// These two events involve exactly the same ball(s).
				// We should keep only the earlier one. Unless the
				// existing event was marked as "maybe" because that
				// means the new event is its re-evaluation.
				if (!existingEvent.isMaybe && existingEvent.lessThanOrEqual(event)) {
					// The event is not earlier than the existing event
					// so we do not add this one.
					return;
				}
				// Remove the existing event with the same ball(s).
				this.remove(existingEvent.index);
			} else if (existingEvent.lessThan(event)) {
				// The event needs to be re-evaluated after the existing
				// event is processed.
				event.isMaybe = true;
			} else {
				// The existing event needs to be re-evaluated after the
				// event is processed.
				existingEvent.isMaybe = true;
			}
		}

		var idx = this.arr.length;
		event.index = idx;
		this.arr.push(event);
		for (var id of event.ballIDs()) {
			var existingEvents = this.ballEvents.get(id);
			if (!existingEvents) {
				existingEvents = new Set();
				this.ballEvents.set(id, existingEvents);
			}
			existingEvents.add(event);
		}

		// This new event is at the bottom of the heap. We may need to bubble
		// it up.
		this.bubbleUp(idx);
		// this.check();
	}

	remove(idx) {
		// Swap this index to the bottom of the array.
		this.swap(idx, this.arr.length-1);

		// Remove the item at the bottom of the array.
		var event = this.arr.pop();
		for (var id of event.ballIDs()) {
			this.ballEvents.get(id).delete(event);
		}

		// If the removed index was at the bottom of the heap then we're done.
		if (idx === this.arr.length) {
			return;
		}

		this.bubbleUp(idx);
		this.bubbleDown(idx);
	}

	removeMin() {
		var minEvent = this.arr[0];
		this.remove(0);
		return minEvent;
	}

	swap(i, j) {
		if (i === j) {
			return;
		}

		var eventI = this.arr[i];
		var eventJ = this.arr[j];
		eventI.index = j;
		eventJ.index = i;
		this.arr[i] = eventJ;
		this.arr[j] = eventI;
	}

	bubbleUp(idx) {
		var parentIdx = parentIndexOf(idx);
		while (idx > 0 && this.arr[idx].lessThan(this.arr[parentIdx])) {
			// This is less than the parent, need to swap.
			this.swap(idx, parentIdx);
			idx = parentIdx;
			parentIdx = parentIndexOf(idx);
		}
	}

	bubbleDown(idx) {
		var leftIdx = leftChildIndexOf(idx);
		var rightIdx = rightChildIndexOf(idx);
		while (rightIdx < this.arr.length) { // Node has 2 children.
			// Determine which child has the lesser value.
			var lesserIdx = (this.arr[leftIdx].lessThan(this.arr[rightIdx])) ? leftIdx : rightIdx;
			if (this.arr[idx].lessThanOrEqual(this.arr[lesserIdx])) {
				return; // No need to continue.
			}

			// This is greater than the lesser child, need swap.
			this.swap(idx, lesserIdx);
			idx = lesserIdx;
			leftIdx = leftChildIndexOf(idx);
			rightIdx = rightChildIndexOf(idx);
		}

		// Check if there is only 1 (left) child.
		if (leftIdx < this.arr.length && this.arr[leftIdx].lessThan(this.arr[idx])) {
			// Left child is less than this, need to swap down.
			this.swap(idx, leftIdx);
		}
		// There are no more entries further down.
		// this.check();
	}
}

/*
- place balls in quadtree
- for each ball:
	- get bounding rectangle (with radius buffer and trajectory stretch)
	- get other balls centered within that rectangle in the quadtree
	  (be sure not to check against itself)
	- for each other ball:
	 	- determine the time until a collision (if one will occur on the
	 	  current trajectory of the balls)
	 	- create an Event for the collision, referencing the time and the
	 	  belligerents.
	- for any collision events, find the earliest one and put it in a
	  priority min-heap
	  	- need a way to efficiently lookup events in the heap by belligerent
	  	  so that there is only ever one event in the heap per belligerent
	  	- if there is already an event in the minheap for either belligerent
	  	  in the heap, keep only the earlier one
	  	- break ties by (time, b1.ID, b2.ID) where b1.ID < b2.ID
	  	- Reasoning: an earlier collision event implies that the later
	  	  event involving the same belligerent will not actually occur (at
	  	  least not at the same time).
- while timePassed < frameDuration:
	- take an event off of the minheap
	- if the heap is empty or the event time is greater than or equal to the
	  frame duration then break
	- move all balls by (event.time - timePassed)
		- would be great if we could do this via the quadtree
	- perform the event collision
	- for the belligerents involved:
		- use the quadtree again to find other balls centered within bounding
		  rectangle
		- create events for any new collisions, add timePassed to event time
		- also check for wall collisions
		- put earliest event in the minheap, again ensuring there is only one
		  event per belligerent
- move all balls by the remaining (frameDuration - timePassed)
*/

class Simulation {
	constructor() {
		// Place balls in quadtree.
		var ballRadius = 5;
		var rectBuffer = 3*ballRadius;
		var maxSize = 800;

		this.ballTree = new QuadTree(0, 0, maxSize, maxSize);
		this.balls = [];

		for (var i = 0; i < 1000; i++) {
			var x = rectBuffer/2 + Math.random() * (maxSize - rectBuffer);
			var y = rectBuffer/2 + Math.random() * (maxSize - rectBuffer);
			var circle = new Circle(x, y, ballRadius);

			var nearbyCircles = this.ballTree.itemsInRange(circle.boundingRect(ballRadius));
			while (nearbyCircles.length > 0) {
				x = rectBuffer/2 + Math.random() * (maxSize - rectBuffer);
				y = rectBuffer/2 + Math.random() * (maxSize - rectBuffer);
				circle = new Circle(x, y, ballRadius);
				nearbyCircles = this.ballTree.itemsInRange(circle.boundingRect(ballRadius));
			}

			var dx = 2*Math.random()-1;
			var dy = 2*Math.random()-1;
			var vec = new Vector(dx, dy);

			this.balls.push(new BouncingBall(this.ballTree, circle, vec));
		}
		// Infect 3 balls.
		for (var i = 0; i < 3; i++) {
			this.balls[i].health = infected;
			this.balls[i].timeOfInfection = Date.now();
		}

		// Create walls.
		this.walls = [
			new BoundaryWall(new Point(0, 0), new Vector(0, 1), new Rectangle(-rectBuffer, -rectBuffer, maxSize+rectBuffer, rectBuffer)), // Top Wall.
			new BoundaryWall(new Point(0, 0), new Vector(1, 0), new Rectangle(-rectBuffer, -ballRadius, rectBuffer, maxSize+rectBuffer)), // Left Wall.
			new BoundaryWall(new Point(maxSize, maxSize), new Vector(0, -1), new Rectangle(-rectBuffer, maxSize - rectBuffer, maxSize+rectBuffer, maxSize+rectBuffer)), // Bottom Wall.
			new BoundaryWall(new Point(maxSize, maxSize), new Vector(-1, 0), new Rectangle(maxSize - rectBuffer, -rectBuffer, maxSize+rectBuffer, maxSize+rectBuffer)), // Right Wall.
		];

		this.canvas = document.getElementById("sim_area");
		this.drawCtx = this.canvas.getContext('2d');

		this.qTreeToggle = document.getElementById("quadtree");
		this.vectorsToggle = document.getElementById("vectors");
		this.socialDistancingControl = document.getElementById("socdist");

		this.ballRadius = ballRadius;
		this.rectBuffer = rectBuffer;
		this.frameNum = 0;
	}

	draw() {
		this.drawCtx.clearRect(0, 0, this.canvas.width, this.canvas.height);
		for (var ball of this.balls) {
			ball.draw(this.drawCtx);
		}
		if (this.vectorsToggle.isToggled) {
			for (var ball of this.balls) {
				this.drawCtx.stroke(ball.vectorPath2D());
			}
		}
		if (this.qTreeToggle.isToggled) {
			for (var path of this.ballTree.path2Ds()) {
				this.drawCtx.stroke(path);
			}
		}
	}

	getSocDist() {
		switch (this.socialDistancingControl.value) {
		case "moderate":
			return 5.1;
		case "extensive":
			return 7.2;
		default:
			return 0;
		}
	}

	doSocialDistancing() {
		var socDist = this.getSocDist();
		if (!socDist) {
			return;
		}

		var socDistSqrt = Math.sqrt(socDist);
		var maxVSqrd = 4 - 0.3*socDist;
		var repulsorRange = this.ballRadius*6;
		for (var ball of this.balls) {
			// If the ball's velocity magnitude squared is greater than the max
			// allowed then the velocity vector is attenuated.
			var vSqrd = ball.vector.scalarProduct(ball.vector);
			if (vSqrd > maxVSqrd) {
				var scaleFactor = 1 - 0.002*(vSqrd/maxVSqrd);
				ball.vector = ball.vector.scale(scaleFactor);
			}

			var nearbyBalls = this.ballTree.itemsInRange(ball.boundingRect(repulsorRange - this.ballRadius));
			for (var nearbyBall of nearbyBalls) {
				if (nearbyBall.id === ball.id) { continue; }

				var posVec = nearbyBall.circle.vectorFrom(ball.circle);
				var dSqrd = posVec.scalarProduct(posVec);

				if (dSqrd > (repulsorRange*repulsorRange)) { continue ; }

				var uVec = posVec.scale(1/Math.sqrt(dSqrd));
				var dVelVec = uVec.scale(0.03*socDist*ball.mass/dSqrd);
				nearbyBall.vector = nearbyBall.vector.add(dVelVec);
			}
		}
	}

	computeNextFrame() {
		this.frameNum++;

		var frameDuration = 0.25;

		// Create an event heap.
		var eventHeap = new CollisionEventsMinHeap();

		// Add collision events between balls.
		for (var ball of this.balls) {
			for (var nearbyBall of this.ballTree.itemsInRange(ball.boundingRect(this.rectBuffer))) {
				if (nearbyBall.id <= ball.id) {
					// if the ball IDs are the same, don't check.
					// if the nearby ball has a lesser ID it would have
					// already been checked against this ball in the outer
					// loop.
					continue;
				}

				var dt = ball.timeUntilBallCollision(nearbyBall);
				if (dt === null) {
					continue; // No collision.
				}

				if (dt > frameDuration) {
					continue; // Event occurs after this frame ends.
				}

				var detectedCollision = new CollisionEvent(dt, ball, nearbyBall);
				// console.log("Detected Collision:", detectedCollision.key(), detectedCollision);
				eventHeap.add(detectedCollision);
			}
		}

		// Add collision events between walls.
		for (var wall of this.walls) {
			for (var nearbyBall of this.ballTree.itemsInRange(wall.region)) {
				var dt = nearbyBall.timeUntilWallCollision(wall);
				if (dt === null) {
					continue; // No wall collision.
				}

				if (dt > frameDuration) {
					continue; // Event occurs after this frame ends.
				}

				var detectedCollision = new CollisionEvent(dt, wall, nearbyBall);
				// console.log("Detected Collision:", detectedCollision.key(), detectedCollision);
				eventHeap.add(detectedCollision);
			}
		}

		var t = 0; // Set time.
		while (!eventHeap.isEmpty()) {
			var event = eventHeap.removeMin();

			if (event.isMaybe && !event.reevaluate(t, frameDuration)) {
				continue;
			}

			// Move all balls to the point of this event.
			if (event.t < t) {
				console.log("past event", this.frameNum, t, event);
				throw new Error("event is in the past");
			}
			if (event.t > t) {
				for (var ball of this.balls) {
					ball.move(event.t-t);
				}
				t = event.t;
			}

			// console.log("Performing Collision:", event.key(), event);
			event.perform();

			// Create any new future events for the balls involved.
			for (var ball of event.balls()) {
				for (var nearbyBall of this.ballTree.itemsInRange(ball.boundingRect(this.rectBuffer))) {
					if (nearbyBall.id === ball.id) {
						continue; // Skip self.
					}

					var dt = ball.timeUntilBallCollision(nearbyBall);
					if (dt === null) {
						continue; // No collision.
					}

					if (t+dt > frameDuration) {
						continue; // Event occurs after this frame ends.
					}

					var detectedCollision = new CollisionEvent(t + dt, ball, nearbyBall);
					// console.log("Detected Another Collision:", detectedCollision.key(), detectedCollision);
					eventHeap.add(detectedCollision);
				}

				// Check for any wall collisions.
				for (var wall of this.walls) {
					var dt = ball.timeUntilWallCollision(wall);
					if (dt === null) {
						continue; // No wall collision.
					}

					if (t + dt > frameDuration) {
						continue; // Event occurs after this frame ends.
					}

					var detectedCollision = new CollisionEvent(t + dt, ball, wall);
					// console.log("Detected Another Collision:", detectedCollision.key(), detectedCollision);
					eventHeap.add(detectedCollision);
				}
			}
		}

		// No more events. Perform a final move of all balls and draw.
		var dt = frameDuration - t;
		if (dt > 0) {
			for (var ball of this.balls) {
				ball.move(dt);
			}
		}

		this.doSocialDistancing();
	}

	run() {
		this.computeNextFrame();
		this.draw();
		this.timeoutID = setTimeout(function(simulation) { simulation.run(); }, 5, this);
	}

	cancel() {
		console.log("cancelling");
		clearTimeout(this.timeoutID);
	}
}

function run() {
	var simulation = new Simulation();
	// console.log(simulation);

	var resetButton = document.getElementById("reset");
	resetButton.addEventListener("click", function() {
		simulation.cancel();
		simulation = new Simulation();
		console.log("running simulation");
		simulation.run();
	});

	console.log("running simulation");
	simulation.run();
};

document.addEventListener('readystatechange', (event) => {
    if (document.readyState === "complete") {
    	console.log("READY");
    	run();
    }
});
